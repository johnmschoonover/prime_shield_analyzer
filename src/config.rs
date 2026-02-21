use clap::Parser;

/// A high-performance Rust tool for analyzing structural bias in consecutive prime sums.
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
pub struct Config {
    /// The upper bound N, expressed as an exponent for 10^N.
    /// E.g., if --max-exponent 10 is provided, N will be 10^10.
    #[arg(short = 'E', long)]
    pub max_exponent: u32,

    /// Number of resolution bins for the time-series output.
    #[arg(short, long, default_value_t = 1000)]
    pub bins: usize,

    /// Directory for output files.
    #[arg(short, long, default_value = "results")]
    pub output_dir: String,
    /// Manually set the sieve segment size in Kilobytes (KB).
    #[arg(long, default_value_t = 128)]
    pub segment_size_kb: usize,

    /// Minimum chunk size, in u64 words, for Rayon parallel sieve chunking.
    #[arg(long, default_value_t = 4096)]
    pub sieve_min_chunk_words: usize,

    /// Multiplier used to derive number of chunks from thread count.
    #[arg(long, default_value_t = 4)]
    pub sieve_chunk_multiplier: usize,

    /// A comma-separated list of prime gap sizes to track (e.g., "2,4,6,12"). All gaps must be even and > 0.
    #[arg(long, default_value = "2,4,6,12,30,34", value_delimiter = ',')]
    pub gaps: Vec<u64>,

    /// Generate a self-contained HTML report with interactive charts.
    #[arg(long)]
    pub web_report: bool,

    /// Run a sieve-only benchmark sweep for tuning and exit.
    #[arg(long)]
    pub benchmark_tuning: bool,

    /// Segment sizes (KB) to benchmark when --benchmark-tuning is enabled.
    #[arg(long, default_value = "1024,4096,8192,16384", value_delimiter = ',')]
    pub benchmark_segment_kb_list: Vec<usize>,

    /// Min-chunk sizes (u64 words) to benchmark when --benchmark-tuning is enabled.
    #[arg(long, default_value = "256,512,1024,4096", value_delimiter = ',')]
    pub benchmark_min_chunk_words_list: Vec<usize>,

    /// Number of measured segments per benchmark configuration.
    #[arg(long, default_value_t = 8)]
    pub benchmark_segments: usize,

    /// Number of warmup segments per benchmark configuration.
    #[arg(long, default_value_t = 2)]
    pub benchmark_warmup_segments: usize,

    /// Disable batch sieve precomputation (ensure_range) for A/B performance testing.
    #[arg(long)]
    pub disable_precompute: bool,

    /// Number of primes per processing batch. Set to 1 for near-streaming behavior.
    #[arg(long, default_value_t = 262_144)]
    pub batch_size: usize,

    /// Disable progress bar rendering for cleaner performance measurements.
    #[arg(long)]
    pub no_progress: bool,
}
