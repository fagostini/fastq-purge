import argparse
import logging
import os
import pathlib
import re
import sys
import tempfile
from collections import defaultdict
from hashlib import sha256
from importlib.metadata import version
from itertools import repeat
from pickle import dumps

import dnaio
import multiprocess as mp
import polars
import psutil
from loky import get_reusable_executor
from multiprocess import Process, Semaphore
from multiprocess.managers import BaseManager
from multiprocess.pool import Pool
from rbloom import Bloom
from rich.logging import RichHandler
from rich.progress import track
from rich_argparse import ArgumentDefaultsRichHelpFormatter

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(message)s",
    datefmt="%H:%M:%S",
    force=True,
    handlers=[
        RichHandler(rich_tracebacks=True, tracebacks_show_locals=True, markup=True)
    ],
)
# Create a local logger for the module
_logger = logging.getLogger(__name__.split(".")[0])

# Define a global variable to store the undetermined set
# This is needed to avoid passing the object to the child process
undetermined_set = set()


def parse_args() -> argparse.Namespace:
    """
    Parse command line arguments.

    Returns:
        argparse.Namespace: Parsed arguments.
    """
    parser = argparse.ArgumentParser(
        prog=f"{sys.argv[0].split('/')[-1]}",
        usage="%(prog)s [options]",
        description="Package to purge fastq files from unwanted reads",
        formatter_class=ArgumentDefaultsRichHelpFormatter,
    )
    parser.add_argument(
        "--flowcell-path",
        type=pathlib.Path,
        help="""Path to the flowcell directory. It is used to extract the lane information from the
        SampleSheet.csv file, and to identify the undetermined and assigned files.
        """,
        default=None,
        required=False,
    )
    parser.add_argument(
        "--undetermined-path",
        type=pathlib.Path,
        help="""Path to the target fastq file(s) to purge. It can be gzipped or not.
        It can be a single file or a directory. If a directory is provided,
        all files in the directory that match the pattern '*.fq*' or '*.fastq*'
        will be used as target files. The search is recursive.""",
        default=None,
        required=False,
    )
    parser.add_argument(
        "--output-path",
        type=pathlib.Path,
        default=None,
        help="""Path to the output folder where the purged fastq files will be saved.
        If not provided, the output files will be saved in the same directory as the undetermined files.
        If the output path does not exist, it will be created.""",
        required=False,
    )
    parser.add_argument(
        "--assigned-path",
        type=pathlib.Path,
        default=None,
        help="""Path to the bloom filter sources file(s). They can be gzipped or not.
        It can be a single file or multiple files separated by spaces, or a directory.
        If a directory is provided, all files in the directory that match the pattern
        '*.fq*' or '*.fastq*' will be used as bloom filter sources. The search is recursive.""",
        required=False,
        nargs="+",
    )
    parser.add_argument(
        "--sample-sheet",
        type=pathlib.Path,
        default=None,
        help="""Path to the sample sheet file. It is used to extract the sample and lane
        information, and to match the undetermined files with the assigned files.
        """,
        required=False,
    )
    parser.add_argument(
        "--keep-original",
        action="store_true",
        help="""Keep the original undetermined fastq files after purging.
        If this option is not set, the original files will be removed after purging,
        while if set, the purged files will be saved with the suffix '.purged' added to the original file basename.
        """,
    )
    parser.add_argument(
        "--method",
        type=str,
        default="exact",
        choices=["exact", "approx"],
        help="""Method to use for filtering. The 'exact' method uses a python set to store 
        the undetermined reads. The 'approx' method uses a bloom filter to store the undetermined
        reads. The 'exact' method is more resource intensive, but it is more accurate. The
        'approx' method is less resource intensive, but it has a false positive rate associated with it.
        """,
    )
    parser.add_argument(
        "--max-items",
        type=int,
        default=1000000000,
        help="""The estimated number of items in the bloom filter. It is used together with the
        false positive rate to calculate the size of the bloom filter. The default is 1 billion.
        This is ignored if the method is 'exact'.""",
    )
    parser.add_argument(
        "--fpr",
        type=float,
        default=0.0001,
        help="False positive rate for the bloom filter. The default is 0.0001.",
    )
    parser.add_argument(
        "--threading-method",
        type=str,
        default="mp_pool",
        choices=["loky", "mp_pool", "mp_manager"],
        help="""Threading library and method to use for processing (used for testing purposes).
        'loky' uses the loky library, which is a robust, cross-platform and cross-version
        implementation of the ProcessPoolExecutor class of concurrent.futures. 'mp_pool' uses
        the multiprocess library, which is a fork of the multiprocessing library with enhanced
        serialization using the dill library. 'mp_manager' uses the multiprocess library with a custom
        manager to share data between processes. The 'mp_pool' and 'mp_manager' methods are
        equivalent, but the 'mp_pool' method less resource intensive.""",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Number of threads to use for processing",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Set the logging level",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {version('fastq-purge')}",
        help="Show version and exit",
    )
    return parser.parse_args()


def parse_sample_sheet(sample_sheet: pathlib.Path) -> dict:
    """
    Parse the sample sheet file and return a dictionary of sample names and their corresponding lanes.

    Args:
        sample_sheet (pathlib.Path): Path to the sample sheet file.

    Returns:
        dict: Dictionary of sample names and their corresponding lanes.
    """
    # Read the sample sheet file and stores all lines in a list
    with open(sample_sheet) as input_file:
        lines = input_file.readlines()

    # Find the line where the FCID starts, where the actual data in CSV format starts
    skip_lines = [i for i, line in enumerate(lines) if line.startswith("FCID")]
    if not skip_lines:
        _logger.error(
            "No valid SampleSheet found! 'FCID' header not found in the file. "
        )
        exit(1)
    else:
        skip_lines = skip_lines[0]

    # Write the lines to a temporary file, skipping the header lines
    tmp = tempfile.NamedTemporaryFile()
    with open(tmp.name, "w") as f:
        for line in lines[skip_lines:]:
            _ = f.write(line)

    # Read the temporary file using polars and process the data
    data = (
        polars.read_csv(tmp.name)
        .with_columns(polars.col("Lane").cast(polars.String).str.zfill(3))
        .with_columns(polars.col("Lane").str.pad_start(4, "L"))
        .with_columns(
            [
                polars.col("index")
                .fill_null("")
                .str.len_chars()
                .alias("index1_length"),
                polars.col("index2")
                .fill_null("")
                .str.len_chars()
                .alias("index2_length"),
            ]
        )
        .with_columns(
            polars.concat_str(
                [
                    polars.col("Recipe"),
                    polars.col("index1_length"),
                    polars.col("index2_length"),
                ],
                separator="-",
            ).alias("Recipe")
        )
    )

    # Select only the lanes that need to be deduplicated
    data = (
        data.select(["Lane", "Sample_Project", "Sample_Name"])
        .group_by(["Lane", "Sample_Project"])
        .all()
        .join(
            data.select(["Lane", "Recipe"])
            .unique()
            .group_by(["Lane"])
            .n_unique()
            .sort("Lane")
            .filter(polars.col("Recipe") > 1),
            on="Lane",
            how="semi",
        )
        .sort(["Lane", "Sample_Project"])
    )

    # Create a dictionary to store the sample sheet data
    sample_sheet_dict = defaultdict(dict)
    for lane in data.get_column("Lane").unique(maintain_order=True):
        sample_sheet_dict.setdefault(lane, defaultdict(list))
        for proj, id in (
            data.filter(polars.col("Lane") == lane)
            .select("Sample_Project", "Sample_Name")
            .iter_rows()
        ):
            sample_sheet_dict[lane].setdefault(proj, []).append(id)

    return sample_sheet_dict


def validate_args(args: argparse.Namespace) -> argparse.Namespace:
    """
    Validate command line arguments.

    Args:
        args (argparse.Namespace): Parsed command line arguments.
    Returns:
        argparse.Namespace: Validated arguments.
    """
    # Check if the flowcell path is provided and is a directory
    if args.flowcell_path and args.flowcell_path.is_dir():
        _logger.debug("Flowcell path provided! Validating paths...")
        args.sample_sheet = (
            args.sample_sheet
            or [
                path
                for path in args.flowcell_path.iterdir()
                if re.match(r"^SampleSheet.(csv|txt)$", path.name)
            ][-1]
        )
        args.undetermined_path = args.undetermined_path or args.flowcell_path
        args.assigned_path = args.assigned_path or args.flowcell_path
    elif not args.flowcell_path.is_dir():
        _logger.error(
            f"Flowcell path '{args.flowcell_path}' is not a valid directory! "
            "Please provide a valid flowcell path."
        )
        exit(1)
    else:
        _logger.info("No flowcell path provided, using manually provided arguments")

    if args.sample_sheet is None or not args.sample_sheet.is_file():
        _logger.error(
            "No sample sheet provided or found! Please provide a valid sample sheet file."
        )
        exit(1)
    elif args.undetermined_path is None or not args.undetermined_path.exists():
        _logger.error(
            "No undetermined path provided or found! Please provide a valid undetermined fastq file."
        )
        exit(1)
    elif args.assigned_path is None or not (
        all([x.exists() for x in args.assigned_path])
        if isinstance(args.assigned_path, list)
        else args.assigned_path.exists()
    ):
        _logger.error(
            "No assigned path provided or found! Please provide assigned fastq file(s) or directory."
        )
        exit(1)
    else:
        _logger.debug("Success! All paths are valid.")
    return args


def process_args(args: argparse.Namespace) -> tuple[argparse.Namespace, dict]:
    """Process and modify the command line arguments.
    Args:
        args (argparse.Namespace): Parsed command line arguments.
    Returns:
        argparse.Namespace: Processed arguments with additional information.
    """

    def explode_path_to_fastq_files(
        path: pathlib.Path,
        recursive: bool = False,
        include_regex: str = None,
        exclude_regex: str = None,
    ) -> list[pathlib.Path]:
        """
        Explode a path into a list of files.
        If the path is a directory, it will return all files in the directory that match
        the pattern '*.fq*' or '*.fastq*'. If the path is a file, it will return the file itself.

        Args:
            path (pathlib.Path): Path to the target file or directory.
            recursive (bool): Whether to search recursively in the directory.
        Returns:
            list[pathlib.Path]: List of files.
        """
        # Patterns for either root-only or recursive search
        patterns = ["**/*.fq*", "**/*.fastq*"] if recursive else ["*.fq*", "*.fastq*"]
        files = []
        if path.is_dir():
            files = [f for pattern in patterns for f in path.glob(pattern)]
        elif path.is_file():
            files.append(path)
        else:
            _logger.error(f"Path '{path}' is neither a file nor a directory")
            exit(1)
        if include_regex:
            files = [f for f in files if re.search(include_regex, f.name)]
        if exclude_regex:
            files = [f for f in files if not re.search(exclude_regex, f.name)]
        return files

    def create_undetermined_paired_dict(
        path: list[pathlib.Path],
    ) -> dict[str, list[pathlib.Path]]:
        """
        Create a dictionary of paired files from the list of paths.
        The keys are the lane numbers and the values are tuples of paired files.

        Args:
            path (list[pathlib.Path]): List of paths to the target files.
        Returns:
            dict[str, list[pathlib.Path]]: Dictionary of paired files.
        """
        # Create the list of all undetermined fastq files, excluding purged files
        path = [
            x
            for x in explode_path_to_fastq_files(
                path,
                recursive=True,
                include_regex="^Undetermined",
                exclude_regex="purged",
            )
        ]
        # Create the patterns for the undetermined files
        patterns = sorted(
            list(
                {tuple([x.parent, re.sub(r"_R[12]_", "_R[12]_", x.name)]) for x in path}
            )
        )
        # Create a dictionary with the lane numbers as keys and the paired files as values
        path = {
            re.search(r"_L[0-9]{3}_", pattern).group().strip("_"): tuple(
                pathlib.Path(parent).glob(f"{pattern}")
            )
            for parent, pattern in patterns
        }
        # Sort the values in the dictionary or insert None if only one file is present
        return {
            key: tuple(sorted(list(values)))
            if len(values) == 2
            else tuple([values[0], None])
            for key, values in path.items()
        }

    # Extract and group the undetermined files by lane
    args.undetermined_path = create_undetermined_paired_dict(args.undetermined_path)
    if not args.undetermined_path:
        _logger.error("No undetermined files found! Please provide valid paths.")
        exit(1)

    _logger.debug("Undetermined files:")
    for lane, targets in args.undetermined_path.items():
        for target in targets:
            _logger.debug(f"    {lane}: '{target}'")

    # Check whether the output path exists, if not, create it
    args.output_path = args.output_path or pathlib.Path(
        args.undetermined_path[list(args.undetermined_path.keys())[0]][0].parent
    )
    if not args.output_path.is_dir():
        _logger.debug("The output path does not exist! It will be created...")
        args.output_path.mkdir(parents=True, exist_ok=True)
    # Create a dictionary to store the output paths
    output_dict = defaultdict(tuple)
    for key, values in args.undetermined_path.items():
        output_list = []
        for input_file in values:
            if input_file is None:
                output_list.append(None)
            else:
                if input_file.suffix == ".gz":
                    fq_suffix = input_file.with_suffix("").suffix
                    output_list.append(
                        args.output_path.joinpath(
                            input_file.with_suffix("")
                            .with_suffix(f".purged{fq_suffix}.gz")
                            .name
                        )
                    )
                else:
                    fq_suffix = input_file.suffix
                    output_list.append(
                        args.output_path.joinpath(
                            input_file.with_suffix(f".purged{fq_suffix}").name
                        )
                    )
        output_dict.setdefault(key, tuple(output_list))
    args.output_path = output_dict

    _logger.debug("Output files:")
    for lane, targets in args.output_path.items():
        for target in targets:
            _logger.debug(f"    {lane}: '{target}'")

    # Create the list of all assigned fastq files, excluding udetermined files
    args.assigned_path = list(
        set(
            [
                bs
                for assigned_file in (
                    args.assigned_path
                    if isinstance(args.assigned_path, list)
                    else [args.assigned_path]
                )
                for bs in explode_path_to_fastq_files(
                    assigned_file, recursive=True, exclude_regex="Undetermined"
                )
            ]
        )
    )

    # Collect the source files, removing unnecessary files (e.g. retain only one of the read pairs)
    sources_set = set()
    clean_assigned_files = []
    for assigned_file in sorted(args.assigned_path):
        source_basename = re.sub(r"_[IR][0-9]_", "_", assigned_file.name)
        if source_basename in sources_set:
            _logger.debug(
                f"Ignoring '{assigned_file.name}' as another file with the same basename "
                f"('{re.sub(r'(_001)?.f(ast)?q(.gz)?', '', source_basename)}') was found in the same path"
            )
            continue
        sources_set.add(source_basename)
        clean_assigned_files.append(assigned_file)
    args.assigned_path = clean_assigned_files

    # Group the assigned files by lane
    assigned_dict = defaultdict(list)
    for assigned_file in args.assigned_path:
        try:
            key = re.search(r"_L[0-9]{3}_", assigned_file.name).group().strip("_")
        except AttributeError:
            _logger.warning(
                f"File '{assigned_file}' does not match the expected pattern '_L[0-9]{{3}}_'"
            )
            continue
        else:
            assigned_dict.setdefault(key, []).append(assigned_file)
    args.assigned_path = assigned_dict

    _logger.debug("Assigned files:")
    for key, value in args.assigned_path.items():
        for val in value:
            _logger.debug(f"    {key}: '{val}'")

    # Parse the sample sheet to get the undetermined-assigned mapping
    selection = parse_sample_sheet(args.sample_sheet)
    results = defaultdict(dict)
    for key, value in selection.items():
        results.setdefault(key, defaultdict())
        results[key]["undetermined"] = args.undetermined_path.get(key, tuple())
        results[key]["output"] = args.output_path.get(key, tuple())
        results[key]["assigned"] = list()
        for proj, ids in value.items():
            results[key]["assigned"].append(
                [
                    x
                    for x in assigned_dict.get(key, [])
                    if any([x.name.startswith(p) for l in ids for p in l])
                ]
            )
        if not results[key]["assigned"] or not results[key]["undetermined"]:
            _logger.warning(
                f"Lane '{key}' should undergo purging, but no undetermined or assigned reads were found! "
                "Skipping this lane."
            )
            del results[key]
            continue

        _logger.info(
            f"Lane '{key}' has {len(results[key]['undetermined'])} undetermined reads and {len(results[key]['assigned'])} assigned reads files"
        )
        _logger.debug(
            f"Lane '{key}' details:\n"
            f"The following files will be purged: {', '.join([x.name for x in results[key]['undetermined']])}\n"
            f"From reads in the following files: {', '.join([f.name for x in results[key]['assigned'] for f in x])}'"
        )

    if not results:
        _logger.error(
            "No lanes found for purging! Please provide valid paths or a sample sheet."
        )
        exit(1)

    if all(
        not results[key]["assigned"] or not results[key]["undetermined"]
        for key in results
    ):
        _logger.error(
            "No assigned or undetermined reads found! Please provide valid paths."
        )
        exit(1)

    return (args, results)


def memory_usage(logger: logging.Logger, pid: int = None) -> None:
    """
    Print the memory usage of the current process.

    Args:
        logger (logging.Logger): Logger to use for logging the memory usage.
        pid (int, optional): Process ID to check memory usage for. If None, uses the current process.
    """
    vmem = psutil.virtual_memory()
    smem = psutil.swap_memory()
    prefix = f"PID {pid} - " if pid else ""
    logger.debug(f"{prefix}Memory used: {vmem.used / 1024**3:.2f} GB ({vmem.percent}%)")
    logger.debug(f"{prefix}Swap used: {smem.used / 1024**3:.2f} GB ({smem.percent}%)")


def build_bloom_filter(
    bloom_sources: list[pathlib.Path], max_items: int, fpr: float
) -> Bloom:
    """
    Build a bloom filter from the given sources.

    Args:
        bloom_sources (list[pathlib.Path]): List of paths to the bloom filter source files.
        max_items (int): Maximum number of items in the bloom filter.
        fpr (float): False positive rate for the bloom filter.

    Returns:
        Bloom: The constructed bloom filter.
    """

    def hash_func(obj):
        """
        Hash function to use for the bloom filter. It uses SHA256 to hash the object and
        returns the first 16 bytes as an integer.
        This is a simple hash function, but it can be replaced with a more complex one if needed.

        Args:
            obj: The object to hash.
        Returns:
            int: The hash value. Hash function to be used for the bloom filter.
        """
        # Use the first 16 bytes of the SHA256 hash as the hash value
        # This is a simple hash function, but it can be replaced with a more complex one if needed
        # The hash value is converted to an integer
        # using the system byte order and signed=True to allow negative values
        h = sha256(dumps(obj)).digest()
        return int.from_bytes(h[:16], sys.byteorder, signed=True)

    # Create a bloom filter with the given parameters
    bf = Bloom(max_items, fpr, hash_func)
    for bloom_source in bloom_sources:
        _logger.debug(f"Adding '{bloom_source.name}' to bloom filter")
        with dnaio.open(bloom_source) as reader:
            for record in track(
                reader,
                description=f"Parsing '{bloom_source.name}'...",
                total=None,
                transient=True,
            ):
                bf.add(record.name)
    return bf


def process_target_file(
    target: pathlib.Path,
    output_path: pathlib.Path,
    bloom_filter: Bloom,
    threads: int = 1,
) -> None:
    """
    Process the target file and remove reads that are in the bloom filter.

    Args:
        target (pathlib.Path): Path to the target fastq file.
        bloom_filter (Bloom): The bloom filter to use for filtering.
    """

    # Read the target fastq file and write the purged reads to the output file
    # Building the index allows access to the read.raw property, which contains the full read name
    # If the full read name is not needed, it is possible drop the index generation
    total = 0
    valid = 0
    if target.suffix == ".gz":
        fq_suffix = target.with_suffix("").suffix
        output_file = target.with_suffix("").with_suffix(f".purged{fq_suffix}.gz")
    else:
        fq_suffix = target.suffix
        output_file = target.with_suffix(f".purged{fq_suffix}")

    # If provided, use the output path to save the purged fastq file
    if output_path:
        output_file = output_path.joinpath(output_file.name)

    _logger.info(f"Output file: '{output_file}'")
    with dnaio.open(
        output_file, mode="w", compression_level=7, open_threads=(threads + 1) // 2
    ) as writer:
        _logger.info("Reading target fastq file...")
        with dnaio.open(target, open_threads=(threads + 1) // 2) as reader:
            for record in reader:
                if record.name not in bloom_filter:
                    writer.write(record)
                    valid += 1
                total += 1
                if total % 1000000 == 0:
                    _logger.info(f"    Processed {total} reads")

    logging.info(f"    Removed {total - valid} reads")
    logging.info(f"    Kept {valid} reads")
    logging.info(f"    Purged fastq file saved to '{output_file}'")


def build_undetermined_set(path: pathlib.Path) -> set:
    """
    Build a set of undetermined reads from the given path.

    Args:
        path (pathlib.Path): Path to the undetermined fastq file.
    Returns:
        set: Set of undetermined reads.
    """
    _logger.debug(f"Adding '{path.name}' to set")
    with dnaio.open(path) as reader:
        for record in track(
            reader,
            description=f"Parsing '{path.name}'...",
            total=None,
            transient=True,
        ):
            undetermined_set.add(record.name.split(" ")[0])
    return undetermined_set


def get_logger(name: str, log_level: str = "INFO") -> logging.Logger:
    """
    Get a logger with the given name and log level. Used in child processes.

    Args:
        name (str): Name of the logger.
        log_level (str): Log level to set for the logger.
    Returns:
        logging.Logger: The logger with the given name and log level.
    """
    # Get the logger for the given name
    logger = logging.getLogger(name)
    # Set the logging level
    logger.setLevel(log_level)
    return logger


def process_target_set(
    path: pathlib.Path,
    log_level: str = "INFO",
    buffer_size: int = 1000000,
) -> set:
    """
    Process the target file and remove reads that are in the undetermined set.

    Args:
        path (pathlib.Path): Path to the target fastq file.
        undetermined_set (set): Set of undetermined reads.
        buffer_size (int): Size of the buffer to use for processing.
    Returns:
        tuple: Tuple containing the process ID, target file name, and set of already assigned reads.
    """
    # Use the global variable to access the undetermined set
    undetermined_set = os.undetermined_set

    # Create a logger for the process
    logger = get_logger("loky", log_level)
    pid = os.getpid()

    logger.info(f"Processing target file '{path.name}' with pid {pid}")
    already_assigned = set()
    tmp_set = set()
    memory_usage(logger, pid)
    with dnaio.open(path) as reader:
        for i, record in enumerate(reader):
            tmp_set.add(record.name.split(" ")[0])
            if i % buffer_size == 0 and i > 0:
                logger.debug(f"    {pid}: Processed {i} reads")
                tmp_set.intersection_update(undetermined_set)
                already_assigned.update(tmp_set)
                tmp_set.clear()
        logger.debug(f"    {pid}: Processed {i} reads")
        tmp_set.intersection_update(undetermined_set)
        already_assigned.update(tmp_set)
    logger.debug(
        f"Process {pid} ('{path.name}') finished yielding {len(already_assigned)} assigned reads"
    )
    return (pid, path.name, already_assigned)


def loky_process_target_set(
    assigned_list: list[pathlib.Path],
    threads: int = 1,
    log_level: str = "INFO",
) -> set:
    """
    Process the target file and remove reads that are in the undetermined set.

    Args:
        assigned_dict (dict): Dictionary of target files to process.
        threads (int): Number of threads to use for processing.
        log_level (str): Log level to set for the logger.
    Returns:
        set: Set of already assigned reads.
    """

    # Hackish trick to pass a global variable to the child process
    # by mutating a global variable from a module (such as the os module)
    def set_value(value):
        """
        Set the value of the global variable in the child process.

        Args:
            value (set): The value to set for the global variable.
        """
        _logger.info(f"Setting value in child process {os.getpid()}")
        os.undetermined_set = value

    _logger.info("Initializing loky executor...")
    executor = get_reusable_executor(
        max_workers=threads,
        timeout=30,
        kill_workers=True,
        initializer=set_value,
        initargs=(undetermined_set,),
    )
    _logger.info("Purging already assigned reads...")

    results = executor.map(
        process_target_set,
        assigned_list,
        repeat(log_level),
    )
    return set([x for res in results for x in res[2]])


class CustomManager(BaseManager):
    """Custom manager to share data between processes."""

    # nothing
    pass


class MultiprocessingCustom:
    """Custom class to be used with the multiprocessing manager."""

    # constructor
    def __init__(self, data):
        # store the data in the instance
        self.undetermined = data
        self.buffer_size = 1000000
        self.assigned = set()
        self.log = defaultdict(int)

    def intersect_and_update(self, subset) -> int:
        """
        Intersect the subset with the undetermined set and update the assigned set.

        Args:
            subset (set): The subset to intersect with the undetermined set.
        Returns:
            int: The number of reads that were already assigned.
        """
        subset.intersection_update(self.undetermined)
        count = len(subset)
        self.assigned.update(subset)
        subset.clear()
        return count

    # perform the main task
    def task(self, path):
        """
        Perform the main task of the custom class.
        Args:
            path (pathlib.Path): Path to the target fastq file.
        Returns:
            tuple: Tuple containing the process ID, target file name, and number of already assigned reads.
        """
        tmp_set = set()
        counts = 0
        with dnaio.open(path) as reader:
            for i, record in enumerate(reader):
                tmp_set.add(record.name.split(" ")[0])
                if i % self.buffer_size == 0 and i > 0:
                    counts += self.intersect_and_update(tmp_set)
        counts += self.intersect_and_update(tmp_set)
        _logger.debug(
            f"Process ('{path.name}') finished yielding {counts} already assigned reads"
        )
        self.log[path.name] = counts
        return (path.name, counts)

    # view the log
    def view_log(self):
        """View the log of the custom class."""
        for key, value in self.log.items():
            _logger.info(f"{key}: {value}")

    # get all stored values
    def get_assigned(self):
        """Get all stored values."""
        return self.assigned

    # remove all stored values
    def clear_assigned(self):
        """Remove all stored values."""
        return self.assigned.clear()


def manager_work(shared_custom, path, semaphore, log_level="INFO"):
    """
    Custom function to be executed in a child process.

    Args:
        shared_custom (MultiprocessingCustom): Shared custom class instance.
        path (pathlib.Path): Path to the target fastq file.
        semaphore (Semaphore): Semaphore to limit the number of concurrent processes.
        log_level (str): Log level to set for the logger.
    """
    # acquire the semaphore
    semaphore.acquire()
    # create a logger for the process
    logger = get_logger("multiprocessing", log_level)
    pid = os.getpid()
    logger.debug(f"Process {pid} started: {path.name} with {mp.get_start_method()}")
    # call the function on the shared custom instance
    name, number = shared_custom.task(path)
    # return the result
    logger.info(f"Process {pid} finished: {name} {number}")
    # release the semaphore
    semaphore.release()


def pool_task(path, log_level="INFO") -> tuple:
    """
    Custom function to be executed in a child process.

    Args:
        path (pathlib.Path): Path to the target fastq file.
        log_level (str): Log level to set for the logger.
    Returns:
        tuple: Tuple containing the process ID, target file name, and set of already assigned reads.
    """
    logger = get_logger("multiprocessing", log_level)
    pid = os.getpid()
    logger.debug(f"Process {pid} started: {path.name} with {mp.get_start_method()}")
    global undetermined_set

    memory_usage(logger, pid)

    tmp_set = set()
    assigned_set = set()
    with dnaio.open(path) as reader:
        for i, record in enumerate(reader):
            tmp_set.add(record.name.split(" ")[0])
            if i % 1000000 == 0 and i > 0:
                tmp_set.intersection_update(undetermined_set)
                assigned_set.update(tmp_set)
                tmp_set.clear()
        tmp_set.intersection_update(undetermined_set)
        assigned_set.update(tmp_set)
    _logger.debug(
        f"Process {pid} ('{path.name}') finished yielding {len(assigned_set)} already assigned reads"
    )
    return (path.name, assigned_set)


def multiprocessing_process_target_set(
    assigned_list: list[pathlib.Path],
    threads: int = 1,
    log_level: str = "INFO",
    method: str = "pool",  # either "manager" or "pool"
) -> None:
    """
    Process the target file and remove reads that are in the undetermined set.

    Args:
        assigned_dict (dict): Dictionary of target files to process.
        threads (int): Number of threads to use for processing.
        log_level (str): Log level to set for the logger.
        method (str): Method to use for multiprocessing. Either "manager" or "pool".
    Returns:
        set: Set of already assigned reads.
    """
    mp.set_start_method("fork", force=True)

    already_assigned = set()
    memory_usage(_logger)

    if method == "manager":
        _logger.info("Using multiprocessing custom manager...")
        # Method 1: Using a custom manager
        CustomManager.register("MultiprocessingCustom", MultiprocessingCustom)
        with CustomManager() as manager:
            # create a shared custom class instance
            shared_custom = manager.MultiprocessingCustom(undetermined_set)
            memory_usage(_logger)
            semaphore = Semaphore(threads)
            _logger.debug(f"Creating {len(assigned_list)} child processes...")
            processes = [
                Process(
                    target=manager_work,
                    args=(shared_custom, target, semaphore, log_level),
                )
                for target in assigned_list
            ]
            _logger.debug(f"Starting {len(processes)} child processes")
            for process in processes:
                process.start()
            _logger.debug(f"Waiting for {len(processes)} child processes to finish")
            for process in processes:
                process.join()
            _logger.debug("Child process finished")
            already_assigned.update(shared_custom.get_assigned())
            shared_custom.clear_assigned()
    else:
        _logger.info("Using multiprocessing pool...")
        # Method 2: Using a process pool
        _logger.info("Purging already assigned reads...")
        _logger.debug(f"Creating {len(assigned_list)} child processes...")
        # create and configure the process pool
        with Pool(processes=threads) as pool:
            # issue tasks to the process pool
            results = pool.map(pool_task, assigned_list)
        for res in results:
            _logger.debug(
                f"Process finished: {res[0]} with {len(res[1])} already assigned reads"
            )

        already_assigned = set([x for res in results for x in res[1]])
    return already_assigned


def write_purged_fastq(
    path_undetermined: pathlib.Path,
    path_purged: pathlib.Path,
    assigned_set: set,
    keep_original: bool = False,
    threads: int = 1,
) -> None:
    """
    Write the purged fastq file.

    Args:
        path_undetermined (pathlib.Path): Path to the undetermined fastq file.
        path_purged (pathlib.Path): Path to the purged fastq file.
        assigned_set (set): Set of already assigned reads.
        threads (int): Number of threads to use for processing.
    """
    if path_undetermined[1] is None or path_purged[1] is None:
        _logger.info(f"Writing purged fastq file '{path_purged[0].name}'")
        try:
            with dnaio.open(
                path_purged[0],
                mode="w",
                compression_level=7,
                open_threads=(threads + 1) // 2,
            ) as writer:
                with dnaio.open(
                    path_undetermined[0], open_threads=(threads + 1) // 2
                ) as reader:
                    for record in reader:
                        if record.name.split(" ")[0] not in assigned_set:
                            writer.write(record)
        except RuntimeError as e:
            _logger.error(
                f"Error writing purged fastq file '{path_purged[0].name}': {e}"
            )
            exit(1)
        else:
            if not keep_original:
                logging.info(
                    f"Replacing undetermined fastq file '{path_undetermined[0].name}' "
                )
                path_purged[0].replace(path_undetermined[0])
    else:
        _logger.info(
            f"Writing purged fastq files '{path_purged[0].name}' and '{path_purged[1].name}'"
        )
        try:
            with dnaio.open(
                path_undetermined[0],
                path_undetermined[1],
                open_threads=(threads + 1) // 2,
            ) as reader:
                with dnaio.open(
                    path_purged[0],
                    path_purged[1],
                    mode="w",
                    compression_level=7,
                    open_threads=(threads + 1) // 2,
                ) as writer:
                    for r1, r2 in reader:
                        if r1.name.split(" ")[0] not in assigned_set:
                            writer.write(r1, r2)
        except RuntimeError as e:
            _logger.error(
                f"Error writing purged fastq files '{path_purged[0].name}' and '{path_purged[1].name}': {e}"
            )
            exit(1)
        else:
            if not keep_original:
                logging.info(
                    f"Replacing undetermined fastq files '{path_undetermined[0].name}' and "
                    f"'{path_undetermined[1].name}' with purged fastq files "
                )
                # Replace the undetermined files with the purged files
                path_purged[0].replace(path_undetermined[0])
                path_purged[1].replace(path_undetermined[1])


def main() -> None:
    """Main function to run the script."""

    args = parse_args()

    # Set the logging level based on the command line argument
    _logger.setLevel(args.log_level)

    # Validate the command line arguments
    args, matching_table = process_args(validate_args(args))

    for lane, data in matching_table.items():
        _logger.info(f"Processing lane '{lane}'")
        files_undetermined = [data["undetermined"]]
        files_output = data["output"]
        files_assigned = [x for l in data["assigned"] for x in l if x]
        _logger.debug(f"    Undetermined files: {files_undetermined}")
        _logger.debug(f"    Output files: {files_output}")
        _logger.debug(f"    Assigned files: {files_assigned}")

        if args.method == "approx":
            _logger.info("Using bloom filter method")
            _logger.info("Building bloom filter...")
            # Create a bloom filter with the given parameters
            bf = build_bloom_filter(files_undetermined, args.max_items, args.fpr)
            # Log the bloom filter parameters
            _logger.debug(f"Bloom filter size: {bf.size_in_bits} bits")
            # _logger.debug(f"Hash functions: {bf.hash_func}")
            _logger.debug(f"Number of items: {bf.approx_items:.1f}")

            # Process the target files and remove reads that are in the bloom filter
            for target in files_assigned:
                _logger.info(f"Processing target file '{target}'")
                process_target_file(target, files_output, bf, args.threads)
        else:
            _logger.info("Using exact method")
            memory_usage(_logger)
            _logger.info("Building python set from undetermined reads...")
            for un_file_1, un_file_2 in files_undetermined:
                undetermined_set = build_undetermined_set(un_file_1)
                undetermined_count = len(undetermined_set)
                set_size_mb = sys.getsizeof(undetermined_set) / 1024**2
                _logger.info(
                    f"Number of undetermined reads: {undetermined_count} ({set_size_mb:.2f} MB)"
                )

                memory_usage(_logger)
                if args.threading_method == "loky":
                    # Loky executor
                    already_assigned = loky_process_target_set(
                        files_assigned,
                        threads=args.threads,
                        log_level=args.log_level,
                    )
                else:
                    # Multiprocess
                    method = "pool" if args.threading_method == "mp_pool" else "manager"
                    already_assigned = multiprocessing_process_target_set(
                        files_assigned,
                        threads=args.threads,
                        log_level=args.log_level,
                        method=method,
                    )

                _logger.info(
                    f"Found {len(already_assigned)} duplicates in the undetermined file."
                )
                if already_assigned:
                    write_purged_fastq(
                        (un_file_1, un_file_2),
                        files_output,
                        already_assigned,
                        args.keep_original,
                        args.threads,
                    )

                else:
                    _logger.info(
                        "Found no duplicates in the undetermined reads. Exiting..."
                    )

    _logger.info("Done!")


if __name__ == "__main__":
    main()
