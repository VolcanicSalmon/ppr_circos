#include "ofMain.h"
#include "ofApp.h"

// Usage: trackplot [--data-dir <path>]
//   <path> is the base folder holding rh/ and dm/ (defaults to the built-in testdata path).
int main(int argc, char** argv) {
    ofSetupOpenGL(1000, 1040, OF_WINDOW);
    ofApp* app = new ofApp();
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if ((a == "--data-dir" || a == "-d") && i + 1 < argc) app->dataDir = argv[++i];
        else if (app->dataDir.empty()) app->dataDir = a;   // bare positional also works
    }
    ofRunApp(app);
    return 0;
}
