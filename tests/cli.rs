use assert_cmd::prelude::*;
use predicates::prelude::*;
use std::process::Command;
use tempfile::tempdir;

#[test]
fn test_cli_smoke_test() -> Result<(), Box<dyn std::error::Error>> {
    // Create a temporary directory for the output
    let temp_dir = tempdir()?;
    // Ensure the output directory is 'report' to match the application's output path
    let output_dir = temp_dir.path().join("report");

    // Prepare the command
    // Use assert_cmd::cargo::cargo_bin! for robust binary invocation
    let mut cmd = Command::new(assert_cmd::cargo::cargo_bin!("prime_shield_analyzer"));
    cmd.arg("--max-exponent")
        .arg("5") // Use a slightly higher exponent for a more representative test
        .arg("--output-dir") // Corrected argument name
        .arg(output_dir.to_str().unwrap())
        .arg("--web-report");

    cmd.assert()
        .success()
        .stdout(predicate::str::contains("Report generated"));

    // Assert that the output files were created in the correct location
    assert!(output_dir.join("index.html").exists()); // Check for index.html
    assert!(output_dir.join("gap_spectrum.csv").exists());
    assert!(output_dir.join("oscillation_series.csv").exists());
    assert!(output_dir.join("global_stats.csv").exists());

    // Clean up the temporary directory
    temp_dir.close()?;

    Ok(())
}
