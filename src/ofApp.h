#pragma once
#include "ofMain.h"
#include "ofxGui.h"
#include "ofxDropdown.h"
#include "BwReader.h"
#include <memory>
#include <string>
#include <vector>
#include <unordered_map>
#include <thread>
#include <atomic>
#include <mutex>

struct Arc {
    std::string id;
    long long   size = 0;
    float       a0 = 0.f, a1 = 0.f;
    ofColor     color = ofColor::white;
};

struct Row {
    std::string chrom;
    long long   site = 0;
    std::string strand, target, trigger;
    std::string allen, leftOfSRNA;                // shown in the table (may be 'NA')
    std::string targetLower, triggerLower;        // precomputed lowercase for fast search
    long long   regionStart = 0, regionEnd = 0;
    float       bscore = 0.f;
};

struct Origin { std::string chrom; long long pos = 0; long long end = 0; std::string strand; };

// A transcript-region feature (rh_ppregion_transcript_regions.bed) for the GW track lane
struct TxFeat { long long start = 0, end = 0; std::string name, strand; };

// A switchable dataset (menu entry).
struct Dataset {
    std::string name;
    std::vector<std::string> bw;        // fwd, rev
    std::string tsv, origins, fasta;
    std::string clustersOnTx;           // optional: cluster->transcript BED (DM)
    std::string transcript, fai;        // genome/: transcript regions BED + .fa.fai
    std::string bam;                    // sorted+indexed sRNA BAM (GW phase/DC)
    std::string jbrowse;                // JBrowse URL template (app appends &loc=chrom:start-end)
};

class ofApp : public ofBaseApp {
public:
    std::string dataDir;                           // base dir for all datasets (--data-dir)
    std::vector<std::string> bwPaths;
    std::string tsvPath, originsPath, fastaPath;   // fasta = trigger sequences by id

    void setup() override;
    void update() override;
    void draw() override;
    void windowResized(int w, int h) override;
    void keyPressed(int key) override;
    void mousePressed(int x, int y, int button) override;
    void mouseDragged(int x, int y, int button) override;
    void mouseReleased(int x, int y, int button) override;
    void mouseScrolled(int x, int y, float sx, float sy) override;

private:
    enum class Mode { Overview, Focus };

    void openReaders();
    void loadChroms();
    void layoutArcs();
    void loadTrack(int binsPerArc);
    std::vector<float> sumBinnedChrom(const std::string& chrom, int nBins);
    std::vector<float> sumBinned(const std::string& chrom, long long s, long long e, int nBins);

    void loadTsv(const std::string& path);
    void loadOrigins(const std::string& path);
    void loadTxRegions(const std::string& path);   // transcript regions BED -> txRegions
    void loadFasta(const std::string& path);   // trigger id -> sequence
    void loadClusterTx(const std::string& path);  // cluster id -> transcript(s)
    void buildDatasets();
    void switchDataset(int i);
    void onDatasetChanged(std::string& name);   // ofxDropdown parameter listener
    void filterTable(const std::string& chrom);
    void showTargetsOfTrigger(const std::string& trig);   // table -> all rows of this trigger
    void backToChromTable();                              // return to the chromosome table
    void applySearch();                                   // filter viewRows by the search box
    void copySelection();                                 // selected rows -> clipboard (TSV)
    void selectRow(int idx);                 // set selRow + gather its on-ring origins

    void enterFocus(int originIdx);          // zoom into one origin locus
    void enterFocusTarget();                 // zoom into the target slice-site region
    bool targetDotHit(int x, int y) const;   // click test on the yellow target dot
    void rebinFocus();

    void drawRing();
    void drawChords();
    void drawArrows();                       // red arrows at the trigger's origins
    void drawFocus();
    void drawTable();

    // GW: JBrowse-style two-pane genome-window view (top = trigger origin, bottom = zoom window)
    void openGW();
    void drawGW();
    void closeGW() { gwOpen = false; gwDragSplit = false; }
    void drawCovPane(ofRectangle rc, const std::string& title, const std::string& chrom,
                     long long s, long long e, long long mark,
                     const std::vector<float>& bins, float mx, ofColor col);
    int  gwSliceHit(const ofRectangle& pane, const std::string& chrom,
                    long long s, long long e, int x, int y) const;  // row nearest a clicked slice line
    void openJBrowse(const std::string& chrom, long long s, long long e);  // 'JB' -> browser at this window
    void requestGwPhase();                                   // mark phase/DC dirty (debounced)
    std::string queryGwPhase(const std::string& chrom, long long s, long long e);  // runs gw_phase.py

    float layoutTop() const { return ofGetHeight() * 0.62f; }
    void  recomputeLayout();
    float angleAt(int x, int y) const;
    int   arcAtAngle(float ang) const;
    long long arcSize(const std::string& chrom) const;
    bool  genomicToPoint(const std::string& chrom, long long pos, ofPoint& out) const;
    bool  originAngle(const Origin& o, float& ang) const;   // overview ring angle
    bool  focusArrowAngle(float& ang) const;                // focused origin in the zoom
    int   arrowAtPixel(int x, int y) const;                 // which origin-arrow was clicked

    std::vector<Dataset> datasets;
    int   curDataset = 0;
    ofxPanel gui;
    ofParameter<std::string> datasetParam;
    std::unique_ptr<ofxDropdown> datasetDropdown;

    std::vector<std::unique_ptr<BwReader>> readers;
    std::vector<Arc>                arcs;
    std::vector<std::vector<float>> trackBins;
    float trackMax = 1.f;

    std::vector<Row> rows;
    std::unordered_map<std::string, std::vector<int>> rowsByChrom;
    std::unordered_map<std::string, std::vector<int>> rowsByTrigger;
    std::unordered_map<std::string, std::vector<Origin>> origins;
    std::unordered_map<std::string, std::vector<TxFeat>> txRegions;   // chrom -> transcript regions
    std::unordered_map<std::string, std::string> triggerSeq;   // trigger id -> sequence
    std::unordered_map<std::string, std::string> clusterTx;    // cluster id -> transcript(s)
    std::unordered_map<std::string, std::string> clusterTxLower;// lowercased, for fast search
    std::vector<int> tableRows;
    std::string activeChrom;
    int   selRow = -1;
    std::vector<Origin> selOrigins;          // on-ring origins of the selected trigger
    int   tableScroll = 0;
    float rowH = 15.f;

    enum class TableMode { ByChrom, ByTrigger };
    TableMode tableMode = TableMode::ByChrom;
    std::string triggerFilter;               // trigger shown in ByTrigger mode
    bool ctxOpen = false;                    // right-click context menu
    int  ctxX = 0, ctxY = 0;
    std::string ctxTrigger;

    std::string searchQuery;                 // target/transcript search box
    bool  searchFocused = false;
    std::vector<int> viewRows;               // tableRows after the search filter (what's drawn)
    int   copyLo = -1, copyHi = -1;          // drag-selected row range (indices into viewRows)
    std::string toast; float toastUntil = 0; // brief "copied" feedback

    ofPoint center;
    float radius = 300.f;
    float trackDepth = 55.f;
    float rotation = 0.f;
    int   selectedArc = -1;

    // focus (zoom on one origin) — entered by clicking its arrow, kept until Esc
    Mode  mode = Mode::Overview;
    Origin focusOrigin;
    bool  focusTarget = false;               // focus is a target (yellow) vs origin (red)
    long long viewStart = 0, viewEnd = 0;
    std::vector<float> focusBins;
    float focusMax = 1.f;
    int   focusResBins = 1000;
    const long long originPad = 3000;        // zoom window = origin.pos +/- this
    const float fa0 = -PI + 0.35f;
    const float fa1 =  PI - 0.35f;
    const float arrowR = 30.f;

    bool  dragging = false;
    int   pressX = 0, pressY = 0;
    float dragStartAngle = 0.f, rotAtDragStart = 0.f;
    long long viewStartAtDrag = 0;

    // GW genome-window overlay state
    bool  gwOpen = false;
    ofRectangle gwBtn;                       // the "GW" button in the focus header (screen coords)
    ofRectangle gwCloseBtn, gwSplitBar;      // set each frame in drawGW()
    ofRectangle gwTopPane, gwBotPane;        // pane rects, for scroll hit-testing
    bool  gwDragSplit = false;
    float gwSplit = 0.5f;                    // top-pane fraction of the content height
    float gwContentTop = 0.f, gwContentH = 1.f;
    std::string gwTopChrom, gwBotChrom;
    long long gwTopS = 0, gwTopE = 0, gwBotS = 0, gwBotE = 0;
    long long gwTopMark = -1, gwBotMark = -1;
    long long gwTopMarkEnd = -1;             // origin read end (for start-end title)
    std::string gwTopStrand;
    std::vector<float> gwTopBins, gwBotBins;
    float gwTopMax = 1.f, gwBotMax = 1.f;
    static const int gwBinN = 800;

    // GW phase/DC (computed live by gw_phase.py off the sorted BAM)
    std::string bamPath, gwScript, gwPython;
    std::string gwTopPhase, gwBotPhase;      // "phase=.. DC=.." (guarded by gwPhaseMtx)
    std::mutex  gwPhaseMtx;
    std::atomic<bool> gwPhaseBusy{false};
    bool  gwPhasePending = false;
    float gwPhaseChangeT = 0.f;
    std::string gwTrigger;                   // trigger id shown in the top-pane title
    std::vector<int> gwRevealed;             // rows whose slice-site label is toggled on (click to toggle)
    std::string jbUrl;                       // current dataset's JBrowse URL template ("" = no JB button)
    ofRectangle gwTopJB, gwBotJB;            // 'JB' buttons per pane (screen rects)
};
