#include <emscripten/val.h>
#include <emscripten/bind.h>
#include <emscripten/emscripten.h>
#include <array>
#include <string>
#include <vector>
#include <string>
#include <iostream>
#include <ranges>
#include <algorithm>
#include <thread>

// copied from https://github.com/emscripten-core/emscripten/issues/11070#issuecomment-717675128
namespace emscripten {
    namespace internal {

        template <typename T, typename Allocator>
        struct BindingType<std::vector<T, Allocator>> {
        using ValBinding = BindingType<val>;
        using WireType = ValBinding::WireType;

        static WireType toWireType(const std::vector<T, Allocator> &vec) {
            return ValBinding::toWireType(val::array(vec));
        }

        static std::vector<T, Allocator> fromWireType(WireType value) {
            return vecFromJSArray<T>(ValBinding::fromWireType(value));
        }
    };

    template <typename T>
    struct TypeID<T,
            typename std::enable_if_t<std::is_same<
                                      typename Canonicalized<T>::type,
                    std::vector<typename Canonicalized<T>::type::value_type,
                    typename Canonicalized<T>::type::allocator_type>>::value>> {
    static constexpr TYPEID get() { return TypeID<val>::get(); }
};

}  // namespace internal
}  // namespace emscripten

class LogisticMap {
public:
    int iterationsToSteadyState;
    int iterationsToShow;
    int extraSymmetricalSamplesPerPixel;
    double logScalingFactor;
    double linearUpperBoundFactor;
    double rLowerBound;
    double rUpperBound;
    double startingValue;
    double xLowerBound;
    double xUpperBound;
    bool antiAliasingEnabled;
    bool logarithmicShadingEnabled; // alternative is linear shading
    std::array<unsigned char, 4> backgroundRGBA;
    std::array<unsigned char, 4> maxLogisticMapRGBA;
    std::array<unsigned char, 4> minLogisticMapRGBA;

    // This variable will be updated by CalculateLogisticMap when it is finished calculating and by parameter changes
    // through the UI. It is read only by CalculateLogisticMap.
    // This variable will also not be locked until the very moment CalculateLogisticMao changes it.
    // That is fine since a change in this variable will initiate a recalculation anyways in CalculateLogisticMap - no
    // reason to continue calculating something that will need to be recalculated anyways!
    bool needsToRecalculate;

    bool currentlyCalculating;
    int rValuesCalculated;

private:
    std::vector<std::vector<double>> frequencies;
    std::vector<double> maxFrequencies;
    double LogisticFunction(double r, double input) {
        return r * input * (1 - input);
    }

public:
    LogisticMap() {
        iterationsToSteadyState = 1000;
        iterationsToShow = 2000;
        extraSymmetricalSamplesPerPixel = 2;
        logScalingFactor = 5;
        linearUpperBoundFactor = 3;
        rLowerBound = 2.75;
        rUpperBound = 3.99;
        startingValue = 0.5;
        xLowerBound = 0;
        xUpperBound = 1;
        antiAliasingEnabled = true;
        logarithmicShadingEnabled = true; // alternative is linear shading
        backgroundRGBA.at(0) = 255;
        backgroundRGBA.at(1) = 255;
        backgroundRGBA.at(2) = 255;
        backgroundRGBA.at(3) = 0;
        maxLogisticMapRGBA.at(0) = 0;
        maxLogisticMapRGBA.at(1) = 0;
        maxLogisticMapRGBA.at(2) = 0;
        maxLogisticMapRGBA.at(3) = 255;
        minLogisticMapRGBA.at(0) = 0;
        minLogisticMapRGBA.at(1) = 0;
        minLogisticMapRGBA.at(2) = 0;
        minLogisticMapRGBA.at(3) = 0;
        needsToRecalculate = false;
        currentlyCalculating = false;
        rValuesCalculated = 0;
    }
    void calculateLogisticMap(int canvasWidth, int canvasHeight) {
        std::cout << "calculating logistic map\n";
        currentlyCalculating = true;
        rValuesCalculated = 0;
        double rPixelStep = (rUpperBound - rLowerBound) / canvasWidth;
        int widthSamplesPerPixel = 2 * extraSymmetricalSamplesPerPixel + 1;
        std::vector<double> columnFrequencies(canvasHeight, 0);
        std::vector<double> maxFrequencies(canvasWidth, 0);
        std::vector <std::vector<double>> frequencies(canvasWidth, columnFrequencies);
        double scalingFactor = 1.0 / (widthSamplesPerPixel * iterationsToShow);
        for (long int i = 0; i < canvasWidth; i++) {
            if (needsToRecalculate) {
                std::cout << "cancelled calculation of logistic map\n";
                break;
            }
            double pixelR = rLowerBound + rPixelStep * i;
            auto &currentFrequencies = frequencies.at(i);
            auto &currentMaxFrequency = maxFrequencies.at(i);
            int previousFrequenciesIndex;
            if (i != 0) {
                previousFrequenciesIndex = i - 1;
            } else {
                previousFrequenciesIndex = i;
            }
            auto &previousFrequencies = frequencies.at(previousFrequenciesIndex);
            auto &previousMaxFrequency = maxFrequencies.at(previousFrequenciesIndex);
            int nextFrequenciesIndex;
            if (i != canvasWidth - 1) {
                nextFrequenciesIndex = i + 1;
            } else {
                nextFrequenciesIndex = i;
            }
            auto &nextFrequencies = frequencies.at(nextFrequenciesIndex);
            auto &nextMaxFrequency = maxFrequencies.at(nextFrequenciesIndex);
            for (int j = -extraSymmetricalSamplesPerPixel; j <= extraSymmetricalSamplesPerPixel; j++) {
                if (i == 0 && j < 0) { continue; } // don't go below the rLowerBound
                if (i == canvasWidth - 1 && j > 0) { continue; } // don't go above the rUpperBound
                double currentR = pixelR + j * rPixelStep / widthSamplesPerPixel;
                double currentValue = startingValue;
                for (int iteration = 0; iteration < iterationsToSteadyState; iteration++) {
                    currentValue = LogisticFunction(currentR, currentValue);
                }
                for (int iteration = 0; iteration < iterationsToShow; iteration++) {
                    currentValue = LogisticFunction(currentR, currentValue);
                    if (currentValue > xLowerBound && currentValue < xUpperBound) {
                        double subpixelHeight =
                                canvasHeight * (currentValue - xLowerBound) / (xUpperBound - xLowerBound);
                        if (antiAliasingEnabled) {
                            int lowerPixelHeight = std::floor(subpixelHeight);
                            int upperPixelHeight = std::ceil(subpixelHeight);
                            // won't check for lowerPixelHeight == upperPixelHeight since that will be very rare
                            // subpixelHeight - lowerPixelHeight and upperPixelHeight - subpixelHeight might seem switched but they're not
                            // TODO: remove lots of repetition by defining a function
                            if (j < 0) {
                                previousFrequencies.at(lowerPixelHeight) +=
                                        (upperPixelHeight - subpixelHeight) * (-j / (double) widthSamplesPerPixel) *
                                        scalingFactor;
                                if (previousMaxFrequency < previousFrequencies.at(lowerPixelHeight)) {
                                    previousMaxFrequency = previousFrequencies.at(lowerPixelHeight);
                                }
                                currentFrequencies.at(lowerPixelHeight) += (upperPixelHeight - subpixelHeight) *
                                                                           ((widthSamplesPerPixel + j) /
                                                                            (double) widthSamplesPerPixel) *
                                                                           scalingFactor;
                                if (currentMaxFrequency < currentFrequencies.at(lowerPixelHeight)) {
                                    currentMaxFrequency = currentFrequencies.at(lowerPixelHeight);
                                }
                                previousFrequencies.at(upperPixelHeight) +=
                                        (subpixelHeight - lowerPixelHeight) * (-j / (double) widthSamplesPerPixel) *
                                        scalingFactor;
                                if (previousMaxFrequency < previousFrequencies.at(upperPixelHeight)) {
                                    previousMaxFrequency = previousFrequencies.at(upperPixelHeight);
                                }
                                currentFrequencies.at(upperPixelHeight) += (subpixelHeight - lowerPixelHeight) *
                                                                           ((widthSamplesPerPixel + j) /
                                                                            (double) widthSamplesPerPixel) *
                                                                           scalingFactor;
                                if (currentMaxFrequency < currentFrequencies.at(upperPixelHeight)) {
                                    currentMaxFrequency = currentFrequencies.at(upperPixelHeight);
                                }
                            } else if (j == 0) {
                                currentFrequencies.at(lowerPixelHeight) +=
                                        (upperPixelHeight - subpixelHeight) * scalingFactor;
                                if (currentMaxFrequency < currentFrequencies.at(lowerPixelHeight)) {
                                    currentMaxFrequency = currentFrequencies.at(lowerPixelHeight);
                                }
                                currentFrequencies.at(upperPixelHeight) +=
                                        (subpixelHeight - lowerPixelHeight) * scalingFactor;
                                if (currentMaxFrequency < currentFrequencies.at(upperPixelHeight)) {
                                    currentMaxFrequency = currentFrequencies.at(upperPixelHeight);
                                }
                            } else {
                                currentFrequencies.at(lowerPixelHeight) +=
                                        (upperPixelHeight - subpixelHeight) * (j / (double) widthSamplesPerPixel) *
                                        scalingFactor;
                                if (currentMaxFrequency < currentFrequencies.at(lowerPixelHeight)) {
                                    currentMaxFrequency = currentFrequencies.at(lowerPixelHeight);
                                }
                                nextFrequencies.at(lowerPixelHeight) += (upperPixelHeight - subpixelHeight) *
                                                                        ((widthSamplesPerPixel - j) /
                                                                         (double) widthSamplesPerPixel) * scalingFactor;
                                if (nextMaxFrequency < nextFrequencies.at(lowerPixelHeight)) {
                                    nextMaxFrequency = nextFrequencies.at(lowerPixelHeight);
                                }
                                currentFrequencies.at(upperPixelHeight) +=
                                        (subpixelHeight - lowerPixelHeight) * (j / (double) widthSamplesPerPixel) *
                                        scalingFactor;
                                if (currentMaxFrequency < currentFrequencies.at(upperPixelHeight)) {
                                    currentMaxFrequency = currentFrequencies.at(upperPixelHeight);
                                }
                                nextFrequencies.at(upperPixelHeight) += (subpixelHeight - lowerPixelHeight) *
                                                                        ((widthSamplesPerPixel - j) /
                                                                         (double) widthSamplesPerPixel) * scalingFactor;
                                if (nextMaxFrequency < nextFrequencies.at(upperPixelHeight)) {
                                    nextMaxFrequency = nextFrequencies.at(upperPixelHeight);
                                }
                            }
                        } else {
                            int pixelHeight = std::round(subpixelHeight);
                            currentFrequencies.at(pixelHeight)++;
                            if (currentMaxFrequency < currentFrequencies.at(pixelHeight)) {
                                currentMaxFrequency = currentFrequencies.at(pixelHeight);
                            }
                        }
                    }
                }
                rValuesCalculated++;
            }
        }
        this->frequencies = frequencies;
        this->maxFrequencies = maxFrequencies;
        if (needsToRecalculate) {
            needsToRecalculate = false;
            calculateLogisticMap(canvasWidth, canvasHeight);
        }
        std::cout << "calculated logistic map\n";
        currentlyCalculating = false;
    }

    void drawLogisticMap(emscripten::val canvas) {
        auto ctx = canvas.call<emscripten::val>("getContext", emscripten::val("2d"));
        auto canvasWidth = canvas["clientWidth"].as<int>();
        auto canvasHeight = canvas["clientHeight"].as<int>();

        std::cout << "drawing logistic map\n";
        int rangeRGBA[4];
        for (int i = 0; i < 4; i++) {
            rangeRGBA[i] = maxLogisticMapRGBA[i] - minLogisticMapRGBA[i];
        }
        std::vector<unsigned char> data(canvasWidth * canvasHeight * 4, 0);
        double cutoffFrequency = std::exp(-logScalingFactor);
        for (long int i = 0; i < canvasWidth; i++) {
            int pixelIndex = (canvasWidth * (canvasHeight - 1) + i) * 4;
            std::vector<double> currentFrequencies = frequencies.at(i);
            double max = maxFrequencies.at(i);
            double columnScalingFactor;
            if (logarithmicShadingEnabled) {
                columnScalingFactor = -std::log(max) / logScalingFactor;
                // columnScalingFactor = std::log(1 / max) / logScalingFactor
            } else {
                columnScalingFactor = 1.0 / max;
                // columnScalingFactor = 1 / max
            }
            for (int j = 0; j < canvasHeight; j++) {
                double currentFrequency = currentFrequencies.at(j);
                if (currentFrequency == 0) {
                    for (int i = 0; i < 4; i++) {
                        data.at(pixelIndex + i) = backgroundRGBA[i];
                    }
                } else {
                    double clampedPixelScalingFactor;
                    if (logarithmicShadingEnabled) {
                        double pixelScalingFactor = std::log(currentFrequency) / logScalingFactor + columnScalingFactor;
                        // pixelScalingFactor = std::log(currentFrequency / max) / logScalingFactor
                        clampedPixelScalingFactor = std::min(1.0, std::max(0.0, pixelScalingFactor + 1));
                        // pixelScalingFactor scales from [-1, 0] for x in [e^-logScalingFactor, 1] and (-infinity, -1] for x in (0, e^-logScalingFactor] barring anti-aliasing inaccuracy
                        // clampedPixelScalingFactor scales from [0, 1] for x in [e^-logScalingFactor, 1]
                    } else {
                        double pixelScalingFactor = currentFrequency * columnScalingFactor;
                        clampedPixelScalingFactor = std::min(1.0, std::max(0.0, pixelScalingFactor *
                                                                                linearUpperBoundFactor));
                        // pixelScalingFactor = currentFrequency / max
                        // pixelScalingFactor scales from [0, linearUpperBoundFactor] for x in [0, 1] barring anti-aliasing inaccuracy
                    }
                    for (int i = 0; i < 4; i++) {
                        data.at(pixelIndex + i) = std::round(
                                clampedPixelScalingFactor * rangeRGBA[i] + minLogisticMapRGBA[i]);
                    }
                }
                pixelIndex -= canvasWidth * 4;
            }
        }
        // TODO: maybe make emscripten directly interpret a std::vector<char> as a Uint8ClampedArray
        auto data2 = emscripten::val(data);
        auto Uint8ClampedArray = emscripten::val::global("Uint8ClampedArray");
        auto data4 = Uint8ClampedArray.new_(data2);
        auto ImageData = emscripten::val::global("ImageData");
        auto data6 = ImageData.new_(data4, canvas["clientWidth"], canvas["clientHeight"]);
        ctx.call<void>("putImageData", data6, emscripten::val(0), emscripten::val(0));
        std::cout << "drew logistic map\n";
    }
};

bool ManipulateLogisticMap(bool resizeLogisticMap, bool calculateLogisticMap, bool renderLogisticMap) {
    static auto logisticMap = LogisticMap();
    auto document = emscripten::val::global("document");
    auto canvas = document.call<emscripten::val>("getElementById", emscripten::val("canvas-logistic-map"));
    auto canvasWidth = canvas["clientWidth"].as<int>();
    auto canvasHeight = canvas["clientHeight"].as<int>();
    auto progress = document.call<emscripten::val>("getElementById", emscripten::val("progress-logistic-map"));
    auto progressbar = document.call<emscripten::val>("getElementById", emscripten::val("progressbar-logistic-map"));
    auto progressLabel = document.call<emscripten::val>("getElementById", emscripten::val("label-progressbar-logistic-map"));
    int numRValues = canvasWidth * (2 * logisticMap.extraSymmetricalSamplesPerPixel + 1) - 2 * logisticMap.extraSymmetricalSamplesPerPixel;
    // we subtract 2 * logisticMap.extraSymmetricalSamplesPerPixel because of the subpixel edges that go off the
    // centers of the edge pixels of the map, which we don't calculate
    std::array<bool, 3> finishValues = {false, false, false};
    if (resizeLogisticMap) {
        progress.set("max", emscripten::val(numRValues));
        finishValues.at(0) = true;
    }
    if (calculateLogisticMap) {
        // logisticMap.calculateLogisticMap(canvasWidth, canvasHeight);
        if (logisticMap.currentlyCalculating) {
            logisticMap.needsToRecalculate = true;
        } else {
            logisticMap.currentlyCalculating = true;
            // technically redundant since calculateLogisticMap will also set currentlyCalculating to true, but starting
            // the thread is slower so we will manually set this to be true in time for the renderLogisticMap logic
            std::thread calculations(&LogisticMap::calculateLogisticMap, &logisticMap, canvasWidth, canvasHeight);
            calculations.detach();
            progressbar["style"].set("visibility", emscripten::val("visible"));
            finishValues.at(1) = true;
        }
    }
    if (renderLogisticMap) {
        int rValuesCalculated = logisticMap.rValuesCalculated;
        if (logisticMap.currentlyCalculating) {
            progress.set("value", emscripten::val(rValuesCalculated));
            std::string rValuesCalculatedString = "Rendering Logistic Map: " + std::to_string(rValuesCalculated) + "/" + std::to_string(numRValues);
            progressLabel.set("innerHTML", emscripten::val(rValuesCalculatedString));
        } else {
            progressbar["style"].set("visibility", emscripten::val("hidden"));
            logisticMap.drawLogisticMap(canvas);
            finishValues.at(2) = true;
        }
    }
    return finishValues.at(0) && resizeLogisticMap || finishValues.at(1) && calculateLogisticMap || finishValues.at(2) && renderLogisticMap;
}

void RenderLogisticMap(double DOMHighResTimeStamp) {
    bool manipulationWasFinished = ManipulateLogisticMap(false, false, true);
    emscripten::val window = emscripten::val::global("window");
    if (!manipulationWasFinished) {
        window.call<void>("requestAnimationFrame", emscripten::val::module_property("RenderLogisticMap"));
    }
}

void RenderPlot(double DOMHighResTimeStamp)
{
    emscripten::val document = emscripten::val::global("document");
    emscripten::val canvas = document.call<emscripten::val>("getElementById", emscripten::val("canvas-plot"));
    emscripten::val ctx = canvas.call<emscripten::val>("getContext", emscripten::val("2d"));
    std::cout << "plot\n";
}

void InitializeCanvas(emscripten::val canvas, emscripten::val index, emscripten::val array)
{
    emscripten::val window = emscripten::val::global("window");
    canvas.set("width", emscripten::val(window["innerWidth"].as<double>() / 2.0 - 1));
    canvas.set("height", emscripten::val(window["innerHeight"].as<double>() - 100 - 1));
    emscripten::val ctx = canvas.call<emscripten::val>("getContext", emscripten::val("2d"));
    ctx.set("textAlign", emscripten::val("center"));
    ctx.set("textBaseline", emscripten::val("middle"));
    ctx.set("font", emscripten::val("20px Arial"));
}

void InitializeCanvases(emscripten::val event)
{
    emscripten::val document = emscripten::val::global("document");
    document.call<emscripten::val>("querySelectorAll", emscripten::val(".canvas")).call<void>("forEach", emscripten::val::module_property("InitializeCanvas"));
    emscripten::val window = emscripten::val::global("window");
    ManipulateLogisticMap(true, true, false);
    window.call<void>("requestAnimationFrame", emscripten::val::module_property("RenderPlot"));
    window.call<void>("requestAnimationFrame", emscripten::val::module_property("RenderLogisticMap"));
}

void InteractWithCanvas(emscripten::val event)
{
    emscripten::val window = emscripten::val::global("window");
    // window.call<void>("requestAnimationFrame", emscripten::val::module_property("RenderPlot"));
    // window.call<void>("requestAnimationFrame", emscripten::val::module_property("RenderLogisticMap"));
}

void InitializeAllSettings()
{
}

int main()
{
    emscripten::val window = emscripten::val::global("window");
    emscripten::val document = emscripten::val::global("document");
    window.call<void>("addEventListener", emscripten::val("resize"), emscripten::val::module_property("InitializeCanvases"));
    document.call<void>("addEventListener", emscripten::val("mousedown"), emscripten::val::module_property("InteractWithCanvas"));
    document.call<void>("addEventListener", emscripten::val("mouseup"), emscripten::val::module_property("InteractWithCanvas"));
    document.call<void>("addEventListener", emscripten::val("mousemove"), emscripten::val::module_property("InteractWithCanvas"));

    // initialize width and height of the canvas
    InitializeCanvases(emscripten::val::null());

    // retrieve all settings from localStorage and set the appropriate boxes to "checked" and put the appropriate data into preview
    InitializeAllSettings();

    return 0;
}

EMSCRIPTEN_BINDINGS(bindings)\
{\
  emscripten::function("InitializeCanvases", InitializeCanvases);\
  emscripten::function("InitializeCanvas", InitializeCanvas);\
  emscripten::function("RenderLogisticMap", RenderLogisticMap);\
  emscripten::function("RenderPlot", RenderPlot);\
};