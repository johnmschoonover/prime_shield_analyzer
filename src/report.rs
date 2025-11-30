use crate::config::Config;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::error::Error;
use std::fs;
use std::path::Path;

#[derive(Serialize, Deserialize, Debug)]
struct GapSpectrumData {
    gap_size: u64,
    success_rate: f64,
    theoretical_boost: f64,
    shield_score: u32,
    shield_primes: String,
}

pub fn generate_report(config: &Config, max_n: u64) -> Result<(), Box<dyn Error>> {
    let output_dir = &config.output_dir;

    // Read oscillation_series.csv dynamically
    let osc_path = Path::new(output_dir).join("oscillation_series.csv");
    let mut osc_reader = csv::Reader::from_path(osc_path)?;
    let mut osc_data: Vec<BTreeMap<String, serde_json::Value>> = Vec::new();
    for result in osc_reader.deserialize() {
        let record: BTreeMap<String, serde_json::Value> = result?;
        osc_data.push(record);
    }
    let osc_json = serde_json::to_string(&osc_data)?;

    // Read gap_spectrum.csv, now including all fields for the new chart
    let gap_path = Path::new(output_dir).join("gap_spectrum.csv");
    let mut gap_reader = csv::Reader::from_path(gap_path)?;
    let mut gap_data: Vec<GapSpectrumData> = Vec::new();
    for result in gap_reader.deserialize() {
        let record: GapSpectrumData = result?;
        if record.success_rate > 0.0 { // Only include gaps with data
             gap_data.push(record);
        }
    }
    let gap_json = serde_json::to_string(&gap_data)?;

    let html_content = format!(
        r#"
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Prime Sum Analysis Report</title>
    <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
    <style>
        body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif; margin: 0; background-color: #f8f9fa; color: #212529; }}
        .container {{ max-width: 1200px; margin: 2rem auto; padding: 2rem; background-color: #fff; border-radius: 8px; box-shadow: 0 4px 6px rgba(0,0,0,0.1); }}
        h1, h2 {{ text-align: center; color: #343a40; }}
        .summary {{ text-align: center; margin-bottom: 2rem; color: #6c757d; }}
        .chart-container {{ margin-top: 2rem; }}
    </style>
</head>
<body>
    <div class="container">
        <h1>Prime Sum Analysis Report</h1>
        <div class="summary">
            <p><strong>Max N:</strong> {max_n} | <strong>Analysis Bins:</strong> {bins}</p>
        </div>
        
        <div class="chart-container">
            <h2>Theory Verification</h2>
            <canvas id="verificationChart"></canvas>
        </div>

        <div class="chart-container">
            <h2>S=p_n+p_(n+1)-1 Primality Ratio Oscillation</h2>
            <canvas id="oscillationChart"></canvas>
        </div>
        
        <div class="chart-container">
            <h2>Gap Success Rate Spectrum (Gaps <= 60)</h2>
            <canvas id="gapChart"></canvas>
        </div>
    </div>

    <script>
        const oscData = {osc_json};
        const gapData = {gap_json};
        const targetGaps = {target_gaps_json};

        // --- Verification Chart (New) ---
        function calculateLinearRegression(data) {{
            const n = data.length;
            if (n === 0) return {{ m: 0, b: 0 }};

            let sumX = 0, sumY = 0, sumXY = 0, sumXX = 0;
            data.forEach(p => {{
                sumX += p.x;
                sumY += p.y;
                sumXY += p.x * p.y;
                sumXX += p.x * p.x;
            }});

            const m = (n * sumXY - sumX * sumY) / (n * sumXX - sumX * sumX);
            const b = (sumY - m * sumX) / n;
            return {{ m, b }};
        }}

        const verificationData = gapData.map(d => ({{
            x: d.theoretical_boost,
            y: d.success_rate,
            gap: d.gap_size,
            score: d.shield_score,
            primes: d.shield_primes
        }}));

        const regression = calculateLinearRegression(verificationData);
        const trendlineData = verificationData.map(p => ({{ x: p.x, y: regression.m * p.x + regression.b }}));

        new Chart(document.getElementById('verificationChart'), {{
            type: 'scatter',
            data: {{
                datasets: [
                    {{
                        label: 'Gaps',
                        data: verificationData,
                        backgroundColor: verificationData.map(p => {{
                            if (p.gap === 4) return 'rgba(255, 99, 132, 1)'; // Red for Gap 4
                            if (p.gap === 34) return 'rgba(54, 162, 235, 1)'; // Blue for Gap 34
                            return 'rgba(0, 0, 0, 0.3)'; // Default
                        }}),
                        pointRadius: verificationData.map(p => (p.gap === 4 || p.gap === 34) ? 7 : 4),
                    }},
                    {{
                        label: 'Trendline',
                        data: trendlineData,
                        type: 'line',
                        borderColor: 'rgba(75, 192, 192, 1)',
                        borderWidth: 2,
                        pointRadius: 0,
                        tension: 0.1
                    }}
                ]
            }},
            options: {{
                plugins: {{
                    tooltip: {{
                        callbacks: {{
                            label: function(context) {{
                                const d = context.raw;
                                return `Gap: ${{d.gap}} | Boost: ${{d.x.toFixed(2)}} | Rate: ${{d.y.toFixed(3)}} | Score: ${{d.score}} | Primes: ${{d.primes || 'none'}}`;
                            }}
                        }}
                    }}
                }},
                scales: {{
                    x: {{ title: {{ display: true, text: 'Theoretical Boost' }} }},
                    y: {{ title: {{ display: true, text: 'Observed Success Rate' }} }}
                }}
            }}
        }});


        // --- Oscillation Chart ---
        const oscillationDatasets = [ {{ label: 'Ratio S_p / p', data: oscData.map(d => d.ratio_s_p), borderColor: 'rgba(75, 192, 192, 1)', tension: 0.1 }} ];
        const colors = [
            'rgba(255, 99, 132, 0.5)', 'rgba(54, 162, 235, 0.5)', 'rgba(255, 206, 86, 0.5)',
            'rgba(75, 192, 192, 0.5)', 'rgba(153, 102, 255, 0.5)', 'rgba(255, 159, 64, 0.5)'
        ];
        let colorIndex = 0;
        targetGaps.forEach(gap => {{
            const gapKey = `gap_${{gap}}_rate`;
            if (oscData.length > 0 && oscData[0][gapKey] !== undefined) {{
                oscillationDatasets.push({{
                    label: `Gap ${{gap}} Rate`,
                    data: oscData.map(d => d[gapKey]),
                    borderColor: colors[colorIndex % colors.length],
                    hidden: true,
                }});
                colorIndex++;
            }}
        }});
        new Chart(document.getElementById('oscillationChart'), {{
            type: 'line',
            data: {{ labels: oscData.map(d => d.bin_start), datasets: oscillationDatasets }},
            options: {{ scales: {{ y: {{ title: {{ display: true, text: 'Ratio' }} }}, x: {{ title: {{ display: true, text: 'N (Bin Start)' }} }} }} }}
        }});

        // --- Gap Spectrum Chart ---
        new Chart(document.getElementById('gapChart'), {{
            type: 'bar',
            data: {{
                labels: gapData.filter(d=>d.gap_size <= 60).map(d => d.gap_size),
                datasets: [{{
                    label: 'Success Rate',
                    data: gapData.filter(d=>d.gap_size <= 60).map(d => d.success_rate),
                    backgroundColor: 'rgba(153, 102, 255, 0.6)'
                }}]
            }},
            options: {{ scales: {{ y: {{ beginAtZero: true, title: {{ display: true, text: 'Success Rate' }} }}, x: {{ title: {{ display: true, text: 'Gap Size' }} }} }} }}
        }});
    </script>
</body>
</html>
"#,
        max_n = max_n,
        bins = config.bins,
        osc_json = osc_json,
        gap_json = gap_json,
        target_gaps_json = serde_json::to_string(&config.gaps)?,
    );

    let report_path = Path::new(output_dir).join("index.html");
    fs::write(report_path, html_content)?;

    Ok(())
}
