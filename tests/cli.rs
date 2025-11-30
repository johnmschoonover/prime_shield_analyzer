use assert_cmd::prelude::*;
use predicates::prelude::*;
use std::process::Command;
use tempfile::tempdir;

#[test]
fn test_cli_smoke_test() -> Result<(), Box<dyn std::error::Error>> {
    // Create a temporary directory for the output
    let temp_dir = tempdir()?;
    let output_dir = temp_dir.path().join("results");

    // Prepare the command
    let mut cmd = Command::cargo_bin(assert_cmd::pkg_name!())?;
    cmd.arg("--max-exponent")
        .arg("4") // Use a small exponent to run quickly
        .arg("--output-dir")
        .arg(output_dir.to_str().unwrap())
        .arg("--web-report");

    // Run the command and assert success
    cmd.assert()
        .success()
        .stdout(predicate::str::contains("Report generated"));

    // Assert that the output files were created
    assert!(output_dir.exists());
    assert!(output_dir.join("gap_spectrum.csv").exists());
    assert!(output_dir.join("report.html").exists());

    // Clean up the temporary directory
    temp_dir.close()?;

    Ok(())
}
