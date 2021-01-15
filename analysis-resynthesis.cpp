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
#include <cmath>      // ::sin()
#include <iostream>
#include <vector>

#include "al/app/al_App.hpp"
#include "al/ui/al_ControlGUI.hpp"
#include "al/ui/al_Parameter.hpp"

// define free functions (functions not associated with any class) here
//

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
  if (b > n)  //
    b = 0;
  // linear interpolation
  return (1 - t) * table.data[a] + t * table.data[b];
}

// a constant global "variable" as an alternative to a pre-processor definition
const double SAMPLE_RATE = 48000.0;
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
    if (phase > 1)  //
      phase -= 1;
    return value;
  }
};

// a class that may be used as a Sine oscillator
//
struct Sine : Phasor {
  double operator()() {  //
    return sine(Phasor::operator()());
  }
};

// suggested entry in a table of data resulting from the analysis of the input
// sound.
struct Entry {
  double frequency, amplitude;
};

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

using namespace al;

struct MyApp : App {
  Parameter background{"background", "", 0.0, "", 0.0f, 1.0f};
  Parameter t{"t", "", 0.0, "", 0.0f, 1.0f};
  ControlGUI gui;

  int N{16};  // the number of sine oscillators to use

  std::vector<Sine> sine;
  std::vector<std::vector<Entry>> data;  // Matrix
  // data[i][j]
  // data[i] // the ith "frame" contains a set of peaks (frequency/amplitude
  // pair)
  //
  // data[a] and data[b] (next frame)
  //       data[t] (a linear interpolation of frames a and b
  //

  MyApp(int argc, char *argv[]) {
    // C++ "constructor" called when MyApp is declared
    //

    // XXX Analysis!
    // - reuse your STFT peaks analysis code here
    // - fill the data structure (i.e., std::vector<std::vector<Entry>> data)
    // - accept and process command line arguments here
    //   + the name of a .wav file
    //   + the number of oscillators N
    // - adapt code from wav-read.cpp
  }

  void onInit() override {
    // called a single time just after the app is started
    //

    sine.resize(N);

    // remove this code later. it's just here to test
    for (int n = 0; n < N; n++) {
      sine[n].frequency(220.0 * (1 + n));
    }
  }

  void onCreate() override {
    // called a single time (in a graphics context) before onAnimate or onDraw
    //

    nav().pos(Vec3d(0, 0, 8));  // Set the camera to view the scene

    gui << background;
    gui << t;
    gui.init();

    // Disable nav control; So default keyboard and mouse control is disabled
    navControl().active(false);
  }

  void onAnimate(double dt) override {
    // called over and over just before onDraw
    t.set(t.get() + dt * 0.03);
    if (t > 1) {
      t.set(t.get() - 1);
    }
  }

  void onDraw(Graphics &g) override {
    // called over and over, once per view, per frame. warning! this may be
    // called more than once per frame. for instance, in the context of 3D
    // stereo viewing, this will be called twice per frame. if 6 screens are
    // attached to this system, then onDraw will be called 6 times per frame.
    //
    g.clear(background);
    //
    //

    // Draw th GUI
    gui.draw(g);
  }

  void onSound(AudioIOData &io) override {
    // called over and over, once per audio block
    //

    // XXX Resynthesis!
    // - add code here
    // - use data from the std::vector<std::vector<Entry>> to adjust the
    // frequency of the N oscillators
    // - use the value of the t parameter to determine which part of the sound
    // to resynthesize
    // - use linear interpolation
    //

    while (io()) {
      float i = io.in(0);  // left/mono channel input (if any);

      // add the next sample from each of the N oscillators
      //
      float f = 0;
      for (int n = 0; n < N; n++) {
        f += sine[n]();  // XXX update this line to effect amplitude
      }
      f /= N;  // reduce the amplitude of the result

      io.out(0) = f;
      io.out(1) = f;
    }
  }

  bool onKeyDown(const Keyboard &k) override {
    int midiNote = asciiToMIDI(k.key());
    return true;
  }

  bool onKeyUp(const Keyboard &k) override {
    int midiNote = asciiToMIDI(k.key());
    return true;
  }
};

int main(int argc, char *argv[]) {
  // MyApp constructor called here, given arguments from the command line
  MyApp app(argc, argv);

  app.configureAudio(48000, 512, 2, 1);

  // Start the AlloLib framework's "app" construct. This blocks until the app is
  // quit (or it crashes).
  app.start();

  return 0;
}
