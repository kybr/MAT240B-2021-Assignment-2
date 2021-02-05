// Assignment 2
//
// - make an Analysis-Resynthesis / Sinusoidal Modeling app
// - place this file (assgnment2.cpp) in a folder in allolib_playground/
// - build with ./run.sh folder/analysis-resynthesis.cpp
// - add your code to this file
//
// -- Karl Yerkes / 2021-01-12 / MAT240B
//
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

#include <algorithm>  // std::sort
#include <algorithm>  // std::sort, std::min
#include <cmath>      // ::sin()
#include <complex>
#include <iostream>
#include <valarray>
#include <vector>

#include "al/app/al_App.hpp"
#include "al/ui/al_ControlGUI.hpp"
#include "al/ui/al_Parameter.hpp"

#define DR_WAV_IMPLEMENTATION
#include "dr_wav.h"

const double PI = 3.141592653589793238460;

typedef std::complex<double> Complex;
typedef std::valarray<Complex> CArray;

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

double dbtoa(double db) { return pow(10.0, db / 20.0); }

// Cooley–Tukey FFT (in-place, divide-and-conquer)
// Higher memory requirements and redundancy although more intuitive
void fft(CArray &x) {
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

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

// a class that encapsulates data approximating one cycle of a sine wave
//
struct SineTable {
  std::vector<double> data;
  SineTable(int n = 16384) {
    data.resize(n);
    for (int i = 0; i < n; i++) {
      data[i] = ::sin(M_PI * 2.0 * i / n);
      // printf("%lf\n", data[i]);
    }
  }
};

// a function that works with the class above to return the value of a sine
// wave, given a value **in terms of normalized phase**. p should be on (0, 1).
//
double sine(double p) {
  static SineTable table;
  int n = table.data.size();
  int a = p * n;
  int b = 1 + a;
  double t = p * n - a;
  if (b >= n)  //
    b = 0;
  // linear interpolation
  return (1 - t) * table.data[a] + t * table.data[b];
}

// a constant global "variable" as an alternative to a pre-processor definition
const double SAMPLE_RATE = 48000.0;
const unsigned BLOCK_SIZE = 512;
//#define SAMPLE_RATE (48000.0)  // pre-processor definition

// a class using the operator/functor pattern for audio synthesis/processing. a
// Phasor or "ramp" wave goes from 0 to 1 in a upward ramping sawtooth shape. it
// may be used as a phase value in other synths.
//
struct Phasor {
  double phase = 0;
  double increment = 0;
  void frequency(double hz) {  //
    increment = hz / SAMPLE_RATE;
  }
  // c++ operator overloading
  // + - * / ()  ... , >> new
  double operator()() {
    double value = phase;
    phase += increment;
    if (phase >= 1.0)  //
      phase -= 1;
    return value;
  }
};

// a class that may be used as a Sine oscillator
//
struct Sine : Phasor {
  double amplitude{0};
  double operator()() {  //
    return sine(Phasor::operator()());
  }
};

// suggested entry in a table of data resulting from the analysis of the input
// sound. renamed Peak; was named Entry.
struct Peak {
  double magnitude, frequency;
};

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

using namespace al;

struct MyApp : App {
  Parameter background{"background", "", 0.0, "", 0.0f, 1.0f};
  Parameter db{"db", "", -60.0, "", -60.0f, 0.0f};
  Parameter t{"t", "", 0.0, "", 0.0f, 1.0f};
  ControlGUI gui;

  int N{16};

  std::vector<Sine> sine;
  std::vector<std::vector<Peak>> data;  // Matrix
  // data[i][j]
  // data[i] // the ith "frame" contains a set of peaks (frequency/amplitude
  // pair)
  //
  // data[a] and data[b] (next frame)
  //       data[t] (a linear interpolation of frames a and b
  //

  MyApp(int argc, char *argv[]) {
    // XXX Analysis!
    // - reuse your STFT peaks analysis code here
    // - fill the data structure (i.e., std::vector<std::vector<Peak>> data)
    // - accept and process command line arguments here
    //   + the name of a .wav file
    //   + the number of oscillators N
    // - adapt code from wav-read.cpp

    if (argc < 2) {
      printf("analysis-resynthesis <.wav> [N oscillators]");
      exit(1);
    }

    std::vector<double> input;

    drwav *pWav = drwav_open_file(argv[1]);
    if (pWav == nullptr) {
      exit(-1);
    }

    float *pSampleData = (float *)malloc((size_t)pWav->totalPCMFrameCount *
                                         pWav->channels * sizeof(float));
    drwav_read_f32(pWav, pWav->totalPCMFrameCount, pSampleData);

    drwav_close(pWav);

    if (pWav->channels == 1)
      for (int i = 0; i < pWav->totalPCMFrameCount; i++)
        input.push_back(pSampleData[i]);  // XXX intermittent crash here !??
    else if (pWav->channels == 2) {
      for (int i = 0; i < pWav->totalPCMFrameCount; i++)
        input.push_back((pSampleData[2 * i] + pSampleData[2 * i + 1]) / 2);
    } else {
      printf("can't handle %d channels\n", pWav->channels);
      exit(1);
    }

    if (argc > 2) {
      N = std::stoi(argv[2]);
    }

    int clipSize = 2048;
    int hopSize = 1024;
    int fftSize = 8192;

    for (int n = 0; n < input.size() - clipSize; n += hopSize) {
      std::vector<double> clip(clipSize, 0.0);
      for (int i = 0; i < clip.size(); i++) {
        //          input      *       Hann window
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

      // XXX this is where we might smooth the spectrum
      // - use a low-pass filter to smooth out the noise
      //

      std::vector<Peak> peak;
      for (int i = 1; i < data.size() / 2; i++) {
        // XXX this is where we might filter out maxima
        // - use a threshold value; is the maxima much greater than the noise
        // floor?
        // - use a "width"; is the maxima the largest within an M-sample window?
        //   + but what should be the value of M?
        //

        // only accept maxima
        //
        if (abs(data[i - 1]) < abs(data[i]))
          if (abs(data[i + 1]) < abs(data[i]))
            peak.push_back({abs(data[i]) / (clip.size() / 2),
                            SAMPLE_RATE * i / data.size()});
      }

      std::sort(peak.begin(), peak.end(), [](Peak const &a, Peak const &b) {
        return a.magnitude > b.magnitude;
      });

      // XXX what to do if the number of peaks is less than N?
      // unlikely for any sophisticated input, but totally possible for
      // pathological cases such as silence or step or a singal click
      //

      peak.resize(N);  // throw away the extras

      if (0) {
        std::sort(peak.begin(), peak.end(), [](Peak const &a, Peak const &b) {
          return a.frequency < b.frequency;
        });
      }

      this->data.emplace_back(peak);
    }

    sine.resize(N);
  }

  void onCreate() override {
    gui << background;
    gui << db;
    gui << t;
    gui.init();
    navControl().active(false);
  }

  void onAnimate(double dt) override {
    //
  }

  void onDraw(Graphics &g) override {
    g.clear(background);
    gui.draw(g);
  }

  void onSound(AudioIOData &io) override {
    // XXX Resynthesis!
    // - add code here
    // - use data from the std::vector<std::vector<Peak>> to adjust the
    // frequency of the N oscillators
    // - use the value of the t parameter to determine which part of the sound
    // to resynthesize
    // - use linear interpolation
    //

    // XXX the code above changes frequency and amplitude immediately based on
    // the "current" value of t. this might lead to sudden jumps in frequency or
    // amplitude with might be audible. instead, we might smooth these jumps,
    // passing these "control" signals through some processor such as Line or
    // OnePole.

    while (io()) {
      // float i = io.in(0);  // left/mono channel input (if any);

      double tmp = t.get();
      tmp += 0.11 / 48000.0;
      if (tmp >= 1.0) {
        tmp -= 1.0;
      }
      t.set(tmp);

      double p = tmp * data.size();
      int a = p;
      int b = 1 + a;
      if (b >= data.size())  //
        b = 0;
      p -= a;

      for (int i = 0; i < N; i++) {
        double frequency =
            (1 - p) * data[a][i].frequency + p * data[b][i].frequency;
        double magnitude =
            (1 - p) * data[a][i].magnitude + p * data[b][i].magnitude;

        sine[i].frequency(frequency);
        sine[i].amplitude = magnitude;
      }

      // add the next sample from each of the N oscillators
      //
      float f = 0;
      for (int n = 0; n < N; n++) {
        f += sine[n]() * sine[n].amplitude;
      }
      // f /= N;  // reduce the amplitude of the result

      f *= dbtoa(db.get());

      io.out(0) = f;
      io.out(1) = f;
    }
  }
};

int main(int argc, char *argv[]) {
  MyApp app(argc, argv);
  app.configureAudio(SAMPLE_RATE, BLOCK_SIZE, 2, 1);
  app.start();
  return 0;
}
