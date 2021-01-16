// Assignment 2
//
// XXX Analysis!
// - find the most prominent N peaks in a sequence of sound clips
// - break the input into a sequence of overlaping sound clips
//   + each sound clip is 2048 samples long
//   + each sound clip should overlap the last by 1024 samples (50% overlap)
//   + use a Hann window to prepare the data before FFT analysis
// - use the FFT to analyse each sound clip
//   + use an FFT size of 8192 (pad the empty data with zeros)
//   + find the frequency and amplitude of each peak in the clip
//     * find maxima, calculate frequency, sort these by amplitude
//   + store the frequency and amplitude of the most prominent N of these
// - test this whole process using several given audio files (e.g.,
// sine-sweep.wav, impulse-sweep.wav, sawtooth-sweep.wav)
//
// -- Karl Yerkes / 2021-01-12 / MAT240B
//

#include <algorithm>  // std::sort, std::min
#include <complex>
#include <iostream>
#include <valarray>
#include <vector>

const double PI = 3.141592653589793238460;
const double SAMPLERATE = 48000.0;

typedef std::complex<double> Complex;
typedef std::valarray<Complex> CArray;

// Cooley–Tukey FFT (in-place, divide-and-conquer)
// Higher memory requirements and redundancy although more intuitive
void fft(CArray& x) {
  const size_t N = x.size();
  if (N <= 1) return;

  // divide
  CArray even = x[std::slice(0, N / 2, 2)];
  CArray odd = x[std::slice(1, N / 2, 2)];

  // conquer
  fft(even);
  fft(odd);

  // combine
  for (size_t k = 0; k < N / 2; ++k) {
    Complex t = std::polar(1.0, -2 * PI * k / N) * odd[k];
    x[k] = even[k] + t;
    x[k + N / 2] = even[k] - t;
  }
}

struct Peak {
  double magnitude, frequency;
};

int main(int argc, char* argv[]) {
  // take N from the command line
  unsigned long N = 16;
  if (argc > 1) {
    N = std::stoi(argv[1]);
  }

  // take the data in
  //
  std::vector<double> input;
  double value;
  while (std::cin >> value) {
    input.push_back(value);
  }

  int clipSize = 2048;
  int hopSize = 1024;
  int fftSize = 8192;

  for (int n = 0; n < input.size() - clipSize; n += hopSize) {
    std::vector<double> clip(clipSize, 0.0);
    for (int i = 0; i < clip.size(); i++) {
      // windowed copy
      clip[i] = input[n + i] * (1 - cos(2 * PI * i / clip.size())) / 2;
    }

    CArray data;
    data.resize(fftSize);
    for (int i = 0; i < clip.size(); i++) {
      data[i] = clip[i];
    }
    for (int i = clip.size(); i < fftSize; i++) {
      data[i] = 0.0;
    }

    fft(data);

    std::vector<Peak> peak;
    for (int i = 1; i < data.size() / 2; i++) {
      if (abs(data[i - 1]) < abs(data[i]))
        if (abs(data[i + 1]) < abs(data[i]))
          peak.push_back(
              {abs(data[i]) / (clip.size() / 2), SAMPLERATE * i / data.size()});
    }

    std::sort(peak.begin(), peak.end(), [](Peak const& a, Peak const& b) {
      return a.magnitude > b.magnitude;
    });

    int M = std::min(N, peak.size()) - 1;
    for (int i = 0; i < M; i++) {
      printf("%lf/%lf, ", peak[i].magnitude, peak[i].frequency);
    }
    printf("%lf/%lf\n", peak[M].magnitude, peak[M].frequency);
  }

  // put your code here!
  //
  // for each sound clip of 2048 samples print information on the the N peaks.
  // print the N frequency/amplitude pairs on a single line. separate the
  // elements of each pair using the slash "/" character and separate each pair
  // using comma "," character. For example:
  //
  // 12.005/0.707,24.01/0.51,48.02/0.3 ...
  // 13.006/0.706,26.01/0.50,52.02/0.29 ...
  //
  // take N as an argument on the command line.
  //

  return 0;
}
