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
#include <utility>
#include <optional>

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
    // TODO: make getters and setters for these to automatically update needsToRecalculate
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
    bool sonificationApplyFourierTransform;
    int rawAudioSampleRateFactor;
    const int maxCalculationStage = 2;
    std::array<unsigned char, 4> backgroundRGBA;
    std::array<unsigned char, 4> maxLogisticMapRGBA;
    std::array<unsigned char, 4> minLogisticMapRGBA;
    int widthSamplesPerPixel;
    int numRValues;

    // This variable will be updated by CalculateLogisticMap when it is finished calculating and by parameter changes
    // through the UI. It is read only by CalculateLogisticMap.
    // This variable will also not be locked until the very moment CalculateLogisticMao changes it.
    // That is fine since a change in this variable will initiate a recalculation anyways in CalculateLogisticMap - no
    // reason to continue calculating something that will need to be recalculated anyways!
    // TODO: add mutexes to these variables to prevent race conditions
    bool needsToRecalculate;

    bool currentlyCalculating;
    int rValuesCalculated;
    int calculationStage;

private:
    int canvasWidth;
    int canvasHeight;
    std::array<double, 2> currentMouseCoordinates;
    std::vector<std::vector<double>> frequencies;
    std::vector<std::vector<double>> dataPoints;
    std::vector<unsigned char> imageData;
    std::vector<std::vector<float>> audioData;
    std::vector<double> maxFrequencies;
    std::optional<emscripten::val> audioCtxWrapper;
    std::vector<std::optional<emscripten::val>> gainNodes;
    // requires std::optional wrappers because these are only permitted to be created after the first user interaction
    bool isCurrentlyPlaying;
    int currentlyPlayingXCoord;

    double LogisticFunction(double r, double input) {
        return r * input * (1 - input);
    }

public:
    void set_current_mouse_coordinates (double x, double y) {
        currentMouseCoordinates.at(0) = x;
        currentMouseCoordinates.at(1) = y;
    }
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
        sonificationApplyFourierTransform = false;
        rawAudioSampleRateFactor = 20;
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
        widthSamplesPerPixel = 2 * extraSymmetricalSamplesPerPixel + 1;
        numRValues = canvasWidth * widthSamplesPerPixel - 2 * extraSymmetricalSamplesPerPixel;
        needsToRecalculate = false;
        currentlyCalculating = false;
        rValuesCalculated = 0;
        calculationStage = 0;
        currentMouseCoordinates.at(0) = 0;
        currentMouseCoordinates.at(1) = 0;
        canvasWidth = 0;
        canvasHeight = 0;
        isCurrentlyPlaying = false;
        currentlyPlayingXCoord = 0;
    }
    void resizeLogisticMap(int newCanvasWidth, int newCanvasHeight) {
        maxFrequencies.resize(newCanvasWidth, 0);
        std::vector<double> columnFrequencies(newCanvasHeight, 0);
        frequencies.resize(newCanvasWidth, columnFrequencies);
        if (newCanvasHeight != canvasHeight) {
            for (auto& currentColumnFrequencyColumn : frequencies) {
                currentColumnFrequencyColumn.resize(newCanvasHeight, 0);
            }
        }
        std::vector<double> columnDataPoints(iterationsToShow, 0);
        dataPoints.resize(newCanvasWidth, columnDataPoints);
        gainNodes.resize(newCanvasWidth);
        canvasWidth = newCanvasWidth;
        canvasHeight = newCanvasHeight;
    }
    void calculateLogisticMap() {
        std::cout << "calculating logistic map\n";
        currentlyCalculating = true;
        rValuesCalculated = 0;
        calculationStage = 1;
        {
            // calculationStage 1: do the actual calculations
            double rPixelStep = (rUpperBound - rLowerBound) / canvasWidth;
            double scalingFactor = 1.0 / (widthSamplesPerPixel * iterationsToShow);
            for (long int i = 0; i < canvasWidth; i++) {
                if (needsToRecalculate) {
                    // TODO: cancelling calculation currently chimerical
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
                auto& currentDataPoints = dataPoints.at(i);
                for (int j = -extraSymmetricalSamplesPerPixel; j <= extraSymmetricalSamplesPerPixel; j++) {
                    if (i == 0 && j < 0) { continue; } // don't go below the rLowerBound
                    if (i == canvasWidth - 1 && j > 0) { continue; } // don't go above the rUpperBound
                    double currentR = pixelR + j * rPixelStep / widthSamplesPerPixel;
                    double currentValue = startingValue;
                    for (int iteration = 0; iteration < iterationsToSteadyState; iteration++) {
                        currentValue = LogisticFunction(currentR, currentValue);
                    }
                    bool measureSonification = (j == 0);
                    for (int iteration = 0; iteration < iterationsToShow; iteration++) {
                        currentValue = LogisticFunction(currentR, currentValue);
                        // sonification
                        if (measureSonification) {
                            currentDataPoints.at(iteration) = currentValue;
                        }
                        // visualization
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
                                                                             (double) widthSamplesPerPixel) *
                                                                            scalingFactor;
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
                                                                             (double) widthSamplesPerPixel) *
                                                                            scalingFactor;
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
            if (needsToRecalculate) {
                needsToRecalculate = false;
                calculateLogisticMap();
            }
        }
        rValuesCalculated = 0;
        calculationStage++;
        {
            // calculationStage 2: convert the calculations into formats ready for sonification and visualization
            int rangeRGBA[4];
            for (int i = 0; i < 4; i++) {
                rangeRGBA[i] = maxLogisticMapRGBA[i] - minLogisticMapRGBA[i];
            }
            std::vector<unsigned char> imageData(canvasWidth * canvasHeight * 4, 0);
            std::vector<float> columnAudioData(iterationsToShow, 0.0);
            std::vector<std::vector<float>> audioData(canvasWidth, columnAudioData);
            double cutoffFrequency = std::exp(-logScalingFactor);
            for (long int i = 0; i < canvasWidth; i++) {
                if (needsToRecalculate) {
                    std::cout << "cancelled calculation of logistic map\n";
                    break;
                }
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
                            imageData.at(pixelIndex + i) = backgroundRGBA[i];
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
                            imageData.at(pixelIndex + i) = std::round(
                                    clampedPixelScalingFactor * rangeRGBA[i] + minLogisticMapRGBA[i]);
                        }
                    }
                    pixelIndex -= canvasWidth * 4;
                }
                if (sonificationApplyFourierTransform) {
                    // TODO
                } else {
                    auto& currentDataPoints = dataPoints.at(i);
                    std::vector<float> currentAudioData(currentDataPoints.begin(), currentDataPoints.end());
                    audioData.at(i) = std::move(currentAudioData);
                }
                rValuesCalculated++;
            }
            if (needsToRecalculate) {
                needsToRecalculate = false;
                calculateLogisticMap();
            } else {
                this->imageData = std::move(imageData);
                this->audioData = std::move(audioData);
            }
        }
        std::cout << "calculated logistic map\n";
        currentlyCalculating = false;
    }

    void drawLogisticMap(emscripten::val canvas) {
        auto ctx = canvas.call<emscripten::val>("getContext", emscripten::val("2d"));
        auto canvasWidth = canvas["clientWidth"].as<int>();
        auto canvasHeight = canvas["clientHeight"].as<int>();

        std::cout << "drawing logistic map\n";
        // TODO: maybe make emscripten directly interpret a std::vector<char> as a Uint8ClampedArray
        auto imageDataJsArray = emscripten::val(imageData);
        auto Uint8ClampedArray = emscripten::val::global("Uint8ClampedArray");
        auto imageDataJsUint8ClampedArray = Uint8ClampedArray.new_(imageDataJsArray);
        auto ImageData = emscripten::val::global("ImageData");
        auto imageDataJsObject = ImageData.new_(imageDataJsUint8ClampedArray, canvas["clientWidth"], canvas["clientHeight"]);
        ctx.call<void>("putImageData", imageDataJsObject, emscripten::val(0), emscripten::val(0));
        std::cout << "drew logistic map\n";
    }

    void sonifyLogisticMap() {
        if (!audioCtxWrapper.has_value()) {
            emscripten::val AudioContext = emscripten::val::global("AudioContext");
            audioCtxWrapper = AudioContext.new_();
        }
        if (isCurrentlyPlaying) {
            auto& currentGainNode = gainNodes.at(currentlyPlayingXCoord).value();
            currentGainNode["gain"].set("value", emscripten::val(0));
        }
        auto& audioCtx = audioCtxWrapper.value();
        int x = currentMouseCoordinates.at(0);
        if (x < 0) {x = 0;}
        if (x > canvasWidth) {x = canvasWidth;}
        isCurrentlyPlaying = true;
        currentlyPlayingXCoord = x;
        auto& currentGainNodeWrapper = gainNodes.at(x);
        if (currentGainNodeWrapper.has_value()) {
            auto& currentGainNode = currentGainNodeWrapper.value();
            currentGainNode["gain"].set("value", emscripten::val(1));
        } else {
            auto audioBuffer = audioCtx.call<emscripten::val>("createBuffer", emscripten::val(1),
                                                              emscripten::val(iterationsToShow *
                                                                              rawAudioSampleRateFactor),
                                                              audioCtx["sampleRate"]);
            emscripten::val channelData = audioBuffer.call<emscripten::val>("getChannelData", emscripten::val(0));
            auto &currentAudioData = audioData.at(x);
            if (sonificationApplyFourierTransform) {
                // TODO: apply Fast Fourier Transform to currentAudioData
            }
            for (int i = 0; i < currentAudioData.size(); i++) {
                auto currentAudioDataPoint = currentAudioData.at(i);
                // effectively divide the sample rate by rawAudioSampleRateFactor so that the oscillations between two values are not outside the range of human hearing
                for (int j = 0; j < rawAudioSampleRateFactor; j++) {
                    channelData.set(std::to_string(rawAudioSampleRateFactor * i + j),
                                    emscripten::val(currentAudioDataPoint));
                }
            }
            auto currentAudioBuffer = audioCtx.call<emscripten::val>("createBufferSource");
            currentAudioBuffer.set("loop", emscripten::val(true));
            currentAudioBuffer.set("buffer", audioBuffer);
            emscripten::val GainNode = emscripten::val::global("GainNode");
            gainNodes.at(x) = GainNode.new_(audioCtx);
            auto& gainNode = gainNodes.at(x).value();
            currentAudioBuffer.call<void>("start");
            currentAudioBuffer.call<void>("connect", gainNode);
            gainNode.call<void>("connect", audioCtx["destination"]);
        }
    }

    void desonifyLogisticMap() {
        if (isCurrentlyPlaying) {
            auto& currentGainNode = gainNodes.at(currentlyPlayingXCoord).value();
            currentGainNode["gain"].set("value", emscripten::val(0));
            isCurrentlyPlaying = false;
        }
    }
};

LogisticMap& GetLogisticMap() {
    static auto logisticMap = LogisticMap();
    return logisticMap;
}

bool ManipulateLogisticMap(bool resizeLogisticMap, bool calculateLogisticMap, bool renderLogisticMap) {
    auto& logisticMap = GetLogisticMap();
    auto document = emscripten::val::global("document");
    auto canvas = document.call<emscripten::val>("getElementById", emscripten::val("canvas-logistic-map"));
    auto progress = document.call<emscripten::val>("getElementById", emscripten::val("progress-logistic-map"));
    auto progressbar = document.call<emscripten::val>("getElementById", emscripten::val("progressbar-logistic-map"));
    auto progressLabel = document.call<emscripten::val>("getElementById", emscripten::val("label-progressbar-logistic-map"));
    int numRValues = logisticMap.numRValues;
    // we subtract 2 * logisticMap.extraSymmetricalSamplesPerPixel because of the subpixel edges that go off the
    // centers of the edge pixels of the map, which we don't calculate
    std::array<bool, 4> finishValues = {false, false, false, false};
    if (resizeLogisticMap) {
        progress.set("max", emscripten::val(numRValues));
        auto canvasWidth = canvas["clientWidth"].as<int>();
        auto canvasHeight = canvas["clientHeight"].as<int>();
        logisticMap.resizeLogisticMap(canvasWidth, canvasHeight);
        finishValues.at(0) = true;
    }
    if (calculateLogisticMap) {
        if (logisticMap.currentlyCalculating) {
            logisticMap.needsToRecalculate = true;
        } else {
            logisticMap.currentlyCalculating = true;
            // technically redundant since calculateLogisticMap will also set currentlyCalculating to true, but starting
            // the thread is slower so we will manually set this to be true in time for the renderLogisticMap logic
            std::thread calculations(&LogisticMap::calculateLogisticMap, &logisticMap);
            calculations.detach();
            progressbar["style"].set("visibility", emscripten::val("visible"));
        }
        finishValues.at(1) = true;
    }
    if (renderLogisticMap) {
        int rValuesCalculated = logisticMap.rValuesCalculated;
        if (logisticMap.currentlyCalculating) {
            progress.set("value", emscripten::val(rValuesCalculated));
            std::string rValuesCalculatedString = "Rendering Logistic Map:\nStage " + std::to_string(logisticMap.calculationStage) + "/" + std::to_string(logisticMap.maxCalculationStage) + "\nrValue " + std::to_string(logisticMap.rValuesCalculated) + "/" + std::to_string(numRValues);
            progressLabel.set("innerHTML", emscripten::val(rValuesCalculatedString));
        } else {
            progressbar["style"].set("visibility", emscripten::val("hidden"));
            logisticMap.drawLogisticMap(canvas);
            finishValues.at(2) = true;
        }
    }
    return finishValues.at(0) && resizeLogisticMap || finishValues.at(1) && calculateLogisticMap || finishValues.at(2) && renderLogisticMap;
}

void RenderLogisticMapOverlay(double DOMHighResTimeStamp) {
    bool manipulationWasFinished = ManipulateLogisticMap(false, false, true);
    emscripten::val window = emscripten::val::global("window");
    if (!manipulationWasFinished) {
        window.call<void>("requestAnimationFrame", emscripten::val::module_property("RenderLogisticMapOverlay"));
    }
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
    emscripten::val boundingBox = canvas["parentElement"].call<emscripten::val>("getBoundingClientRect");
    // This bounding box has parameters that are not limited to being an integer, as opposed to directly calling the width and height of the parent div
    // If the div is not an integer, the extra subpixel should be covered by the border of the div
    // subtract 2 to account for the border of at least 1 px
    canvas.set("width", emscripten::val(std::ceil(boundingBox["width"].as<double>() - 2)));
    canvas.set("height", emscripten::val(std::ceil(boundingBox["height"].as<double>() - 2)));
    emscripten::val ctx = canvas.call<emscripten::val>("getContext", emscripten::val("2d"));
    ctx.set("textAlign", emscripten::val("center"));
    ctx.set("textBaseline", emscripten::val("middle"));
    ctx.set("font", emscripten::val("20px Arial"));
}

void InteractWithLogisticMapCanvas(emscripten::val event)
{
    static bool mouseIsDown;
    std::string eventName = event["type"].as<std::string>();
    if (eventName == "mouseup") {
        // mouseup anywhere on the document
        mouseIsDown = false;
        auto& logisticMap = GetLogisticMap();
        if (!logisticMap.currentlyCalculating) {
            logisticMap.desonifyLogisticMap();
        }
    } else if (eventName == "mousemove" || eventName == "mousedown") {
        // mousemove only on the canvas
        // mousedown anywhere on the document
        if (eventName == "mousedown") {
            mouseIsDown = true;
        }
        // TODO: find a way to make canvas-logistic-map-overlay not hardcoded
        if (eventName == "mousemove" || eventName == "mousedown" && event["target"]["id"].as<std::string>() == "canvas-logistic-map-overlay") {
            auto &logisticMap = GetLogisticMap();
            emscripten::val window = emscripten::val::global("window");
            double offsetX = event["offsetX"].as<double>();
            double offsetY = event["offsetY"].as<double>();
            logisticMap.set_current_mouse_coordinates(offsetX, offsetY);
            if (!logisticMap.currentlyCalculating && mouseIsDown) {
                logisticMap.sonifyLogisticMap();
            }
            // window.call<void>("requestAnimationFrame", emscripten::val::module_property("RenderPlot"));
            // window.call<void>("requestAnimationFrame", emscripten::val::module_property("RenderLogisticMapOverlay"));
        }
    }
}

void InitializeCanvases(emscripten::val event)
{
    emscripten::val document = emscripten::val::global("document");
    document.call<emscripten::val>("querySelectorAll", emscripten::val(".canvas")).call<void>("forEach", emscripten::val::module_property("InitializeCanvas"));
    emscripten::val window = emscripten::val::global("window");
    ManipulateLogisticMap(true, true, false);
    window.call<void>("requestAnimationFrame", emscripten::val::module_property("RenderPlot"));
    window.call<void>("requestAnimationFrame", emscripten::val::module_property("RenderLogisticMap"));
    // TODO: make canvas-logistic-map-overlay not hardcoded
    auto logisticMapCanvasOverlay = document.call<emscripten::val>("getElementById", emscripten::val("canvas-logistic-map-overlay"));
    document.call<void>("addEventListener", emscripten::val("mousedown"), emscripten::val::module_property("InteractWithLogisticMapCanvas"));
    document.call<void>("addEventListener", emscripten::val("mouseup"), emscripten::val::module_property("InteractWithLogisticMapCanvas"));
    logisticMapCanvasOverlay.call<void>("addEventListener", emscripten::val("mousemove"), emscripten::val::module_property("InteractWithLogisticMapCanvas"));
}

void InitializeAllSettings()
{
}

int main()
{
    emscripten::val window = emscripten::val::global("window");
    emscripten::val document = emscripten::val::global("document");
    window.call<void>("addEventListener", emscripten::val("resize"), emscripten::val::module_property("InitializeCanvases"));
    /*
    document.call<void>("addEventListener", emscripten::val("mousedown"), emscripten::val::module_property("InteractWithCanvas"));
    document.call<void>("addEventListener", emscripten::val("mouseup"), emscripten::val::module_property("InteractWithCanvas"));
    document.call<void>("addEventListener", emscripten::val("mousemove"), emscripten::val::module_property("InteractWithCanvas"));
    */
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
  emscripten::function("InteractWithLogisticMapCanvas", InteractWithLogisticMapCanvas);\
};