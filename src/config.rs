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

    /// A comma-separated list of prime gap sizes to track (e.g., "2,4,6,12"). All gaps must be even and > 0.
    #[arg(long, default_value = "2,4,6,12,30", value_delimiter = ',')]
    pub gaps: Vec<u64>,

    /// Generate a self-contained HTML report with interactive charts.
    #[arg(long)]
    pub web_report: bool,
}
