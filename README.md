# Adaptive Signal Processing

This repository contains MATLAB implementations of two adaptive signal processing systems developed for an academic project in ESE 5310:

- **Adaptive Notch Filter** (Part A)
- **Adaptive Equalizer** (Part B)

## Project Overview

### Part A: Adaptive Notch Filter
The adaptive notch filter is designed to attenuate sinusoidal interference from a desired signal using an IIR filter structure decomposed into an FIR filter followed by an all-pole filter. The filter coefficients are updated using the LMS (Least Mean Squares) algorithm.

Key features include:
- Single interfering sinusoid removal
- Handling linearly varying interference
- Two-notch cascade filtering for multiple interference removal

**File:** `adaptive_notch_filter.m`

### Part B: Adaptive Equalizer
The adaptive equalizer compensates for channel distortion and noise using an LMS-based adaptive filter. A known training sequence is used to train the equalizer, which is then tested on a separate dataset to evaluate error rates and decision accuracy.

Key features include:
- LMS-based equalizer coefficient updates
- Analysis of filter order, step size, and SNR impact
- Decision error rate calculations

**File:** `adaptive_equalizer.m`

## Repository Structure
- `adaptive_notch_filter.m`: MATLAB code for Part A - Adaptive Notch Filter
- `adaptive_equalizer.m`: MATLAB code for Part B - Adaptive Equalizer
- `adaptive_dsp_report.pdf`: Project report detailing the methodology, results, and analysis

## Getting Started
### Prerequisites
- MATLAB R2021a or later

### Running the Code
1. Clone this repository:
   ```bash
   git clone [<repo-url>](https://github.com/agowd/adaptive_signal_processing.git)
   ```
2. Open MATLAB and navigate to the repository folder.
3. Run the scripts:
   - For Part A: `adaptive_notch_filter.m`
   - For Part B: `adaptive_equalizer.m`

### Output
- The scripts generate time-domain plots, frequency responses, error convergence graphs, and decision error rates.

## Authors
- Aditya Gowd @agowd
- Varun Trivedi @vtriv

## Acknowledgments
- Dr. Khanna, ESE 5310 Instructor

---

Feel free to reach out for any questions or collaboration opportunities!

