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
#include <stdexcept>
#include <mutex>
#include <numbers>
#include <random>
#include <condition_variable>
#include "fftw-3.3.10/api/fftw3.h"

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

    enum class DefaultEquationType { LOGISTIC, EXPONENTIAL, GAMMA };

    // TODO: make getters and setters for these to automatically update needsToRecalculate
    int functionID;
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
    bool sonificationApplyInverseFourierTransform;
    bool sonificationLogarithmicSampling;
    int rawAudioSampleRateFactor;
    int ifftMinFrequency;
    int ifftMaxFrequency;
    int ifftAudioSamples;
    bool doIfftAudioPeakNormalization;
    double changeChannelAudioLerpTime;
    // ifftAudioSamples cannot be less than twice ifftMaxFrequency
    const int maxCalculationStage = 2;
    std::array<unsigned char, 4> backgroundRGBA;
    std::array<unsigned char, 4> maxLogisticMapRGBA;
    std::array<unsigned char, 4> minLogisticMapRGBA;
    int widthSamplesPerPixel;
    int numRValues;
    std::array<double, 2> currentMouseCoordinates;
    int canvasWidth;
    int canvasHeight;
    double overlayTextMargin;
    bool enableOverlay;
    bool doShowOverlay;
    int currentXCoord;
    std::mutex parameterMutex;
    std::condition_variable parameterConditionVariable;
    std::random_device randomDevice;
    std::optional<int> rngSeed;
    bool doStartGraphAtIteration0;

    // This variable will be updated by CalculateLogisticMap when it is finished calculating and by parameter changes
    // through the UI. It is read only by CalculateLogisticMap.
    // This variable will also not be locked until the very moment CalculateLogisticMao changes it.
    // That is fine since a change in this variable will initiate a recalculation anyways in CalculateLogisticMap - no
    // reason to continue calculating something that will need to be recalculated anyways!
    // TODO: add mutexes to these variables to prevent race conditions
    bool safeToRecalculate;
    bool needsToRecalculate;

    bool currentlyCalculating;
    int iterationsCalculated;
    int numIterations;
    int calculationStage;
    bool isCurrentlyPlaying;

    std::vector<std::vector<double>> frequencies;
    std::vector<std::vector<double>> dataPoints;
    std::vector<unsigned char> imageData;
    // std::vector<double> fftwIO;
    fftw_complex* fftwIn;
    std::vector<double> fftwOut;
    std::vector<std::vector<double>> plotData;
    std::vector<std::vector<double>> audioData;
    std::vector<std::vector<double>> graphData;
    std::vector<double> maxFrequencies;
    std::vector<double> minY;
    std::vector<double> maxY;
    std::optional<emscripten::val> audioCtxWrapper;
    std::vector<std::optional<emscripten::val>> gainNodes;
    // requires std::optional wrappers because these are only permitted to be created after the first user interaction
    fftw_plan fftwPlan;

public:
    void set_iterations_to_steady_state(int iterations) {
        std::unique_lock<std::mutex> lock(parameterMutex);
        parameterConditionVariable.wait(lock, [this]{ return this->safeToRecalculate; });
        iterationsToSteadyState = iterations;
        needsToRecalculate = true;
    }
    void set_iterations_to_show(int iterations) {
        std::unique_lock<std::mutex> lock(parameterMutex);
        parameterConditionVariable.wait(lock, [this]{ return this->safeToRecalculate; });
        iterationsToShow = iterations;
        needsToRecalculate = true;
    }
    void set_extra_symmetrical_samples_per_pixel(int extraSamples) {
        std::unique_lock<std::mutex> lock(parameterMutex);
        parameterConditionVariable.wait(lock, [this]{ return this->safeToRecalculate; });
        extraSymmetricalSamplesPerPixel = extraSamples;
        needsToRecalculate = true;
    }

private:

    double LogisticFunction(double r, double input) {
        switch (functionID) {
            case 0:
                return r * input * (1 - input);
            case 1:
                return r * input * std::exp(-input);
            case 2:
                return r * 1 / std::tgamma(input);
            default:
                return r * input * (1 - input);
        }
    }

public:
    void set_current_mouse_coordinates (double x, double y) {
        currentMouseCoordinates.at(0) = x;
        currentMouseCoordinates.at(1) = y;
    }
    LogisticMap() {
        functionID = 0;
        iterationsToSteadyState = 1000;
        iterationsToShow = 2000;
        extraSymmetricalSamplesPerPixel = 2;
        logScalingFactor = 5;
        linearUpperBoundFactor = 3;
        rLowerBound = 2.75;
        rUpperBound = 4;
        startingValue = 0.5;
        xLowerBound = 0;
        xUpperBound = 1;
        antiAliasingEnabled = true;
        logarithmicShadingEnabled = true; // alternative is linear shading
        sonificationApplyInverseFourierTransform = true;
        sonificationLogarithmicSampling = true; // TODO
        rawAudioSampleRateFactor = 20;
        ifftMinFrequency = 100;
        ifftMaxFrequency = 1000;
        ifftAudioSamples = std::pow(2, 14);
        doIfftAudioPeakNormalization = true;
        changeChannelAudioLerpTime = 0.005;
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
        canvasWidth = 0;
        canvasHeight = 0;
        overlayTextMargin = 15;
        enableOverlay = true;
        doShowOverlay = false;
        currentXCoord = 0;
        needsToRecalculate = false;
        currentlyCalculating = false;
        iterationsCalculated = 0;
        calculationStage = 0;
        currentMouseCoordinates.at(0) = 0;
        currentMouseCoordinates.at(1) = 0;
        isCurrentlyPlaying = false;
        rngSeed = std::nullopt;
    }
    template <typename T> void resizeVector(std::vector<T>& vector, int newSize, T defaultValue) {
        vector.resize(newSize, defaultValue);
        std::fill(vector.begin(), vector.begin() + std::min((int) vector.size(), newSize), defaultValue);
    }
    void resizeLogisticMap(int newCanvasWidth, int newCanvasHeight) {
        if (newCanvasWidth <= 0 || newCanvasHeight <= 0) {
            throw std::invalid_argument("resizeLogisticMap only accepts positive integers");
        }
        resizeVector<double>(maxFrequencies, newCanvasWidth, 0.0);
        resizeVector<double>(minY, newCanvasWidth, newCanvasHeight);
        resizeVector<double>(maxY, newCanvasWidth, 0.0);
        std::vector<double> columnFrequencies(newCanvasHeight, 0.0);
        resizeVector<std::vector<double>>(frequencies, newCanvasWidth, columnFrequencies);
        std::vector<double> columnDataPoints(iterationsToShow, 0.0);
        resizeVector<std::vector<double>>(dataPoints, newCanvasWidth, columnDataPoints);
        resizeVector<unsigned char>(imageData, newCanvasWidth * newCanvasHeight * 4, 0);
        std::vector<double> columnPlotData(newCanvasHeight, 0.0);
        resizeVector<std::vector<double>>(plotData, newCanvasWidth, columnPlotData);
        std::vector<double> columnAudioData;
        if (sonificationApplyInverseFourierTransform) {
            resizeVector<double>(columnAudioData, ifftAudioSamples, 0.0);
            // resizeVector<double>(fftwIO, 2 * ifftAudioSamples, 0.0);
            fftw_free(fftwIn);
            fftwIn = (fftw_complex*) fftw_alloc_complex(2 * ifftAudioSamples);
            resizeVector<double>(fftwOut, 2 * ifftAudioSamples, 0.0);
        } else {
            resizeVector<double>(columnAudioData, iterationsToShow * rawAudioSampleRateFactor, 0.0);
        }
        resizeVector<std::vector<double>>(audioData, newCanvasWidth, columnAudioData);
        gainNodes.resize(newCanvasWidth);
        std::vector<double> columnGraphData;
        if (doStartGraphAtIteration0) {
            resizeVector<double>(columnGraphData, iterationsToSteadyState + iterationsToShow, 0.0);
        } else {
            resizeVector<double>(columnGraphData, iterationsToShow, 0.0);
        }
        resizeVector<std::vector<double>>(graphData, newCanvasWidth, columnGraphData);
        //std::fill(gainNodes.begin(), gainNodes.begin() + std::min(gainNodes.size(), newCanvasWidth), defaultValue);
        // auto newFftwPlan = fftw_plan_r2r_1d(fftwIO.size(), &fftwIO[0], &fftwIO[0], FFTW_HC2R, FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
        auto newFftwPlan = fftw_plan_dft_c2r_1d(2 * ifftAudioSamples, fftwIn, &fftwOut[0], FFTW_ESTIMATE);
        fftw_destroy_plan(fftwPlan);
        fftwPlan = newFftwPlan;
        canvasWidth = newCanvasWidth;
        canvasHeight = newCanvasHeight;
        numRValues = canvasWidth * widthSamplesPerPixel - 2 * extraSymmetricalSamplesPerPixel;
    }
private:
    bool tryCalculateLogisticMap() {
        std::cout << "calculating logistic map\n";
        currentlyCalculating = true;
        iterationsCalculated = 0;
        numIterations = numRValues;
        calculationStage = 1;
        std::mt19937 rng(rngSeed.value_or(randomDevice()));
        std::uniform_real_distribution<> randomIfftPhase(0, 2 * std::numbers::pi);
        {
            // calculationStage 1: do the actual calculations
            double rPixelStep = (rUpperBound - rLowerBound) / canvasWidth;
            double scalingFactor = 1.0 / (widthSamplesPerPixel * iterationsToShow);
            for (int i = 0; i < canvasWidth; i++) {
                if (needsToRecalculate) {
                    // TODO: cancelling calculation currently chimerical
                    break;
                }
                double pixelR = rLowerBound + rPixelStep * i;
                {
                    std::scoped_lock lock(parameterMutex);
                    safeToRecalculate = false;
                    auto &currentFrequencies = frequencies.at(i);
                    auto &currentPlotData = plotData.at(i);
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
                    auto &currentDataPoints = dataPoints.at(i);
                    auto &currentAudioData = audioData.at(i);
                    auto &currentGraphData = graphData.at(i);
                    double audioDataIncrement = 0.5 / iterationsToShow;
                    // this is half of the inverse of iterationsToShow since the IFFT is doubled
                    for (int j = -extraSymmetricalSamplesPerPixel; j <= extraSymmetricalSamplesPerPixel; j++) {
                        if (i == 0 && j < 0) { continue; } // don't go below the rLowerBound
                        if (i == canvasWidth - 1 && j > 0) { continue; } // don't go above the rUpperBound
                        double currentR = pixelR + j * rPixelStep / widthSamplesPerPixel;
                        double currentValue = startingValue;
                        for (int iteration = 0; iteration < iterationsToSteadyState; iteration++) {
                            currentValue = LogisticFunction(currentR, currentValue);
                            // normalize currentValue's range from [xLowerBound, xUpperBound) to [0, 1)
                            double normalizedCurrentValue =
                                    (currentValue - xLowerBound) / (xUpperBound - xLowerBound);
                            // graph
                            if (doStartGraphAtIteration0) {
                                currentGraphData.at(iteration) = normalizedCurrentValue;
                            }
                        }
                        bool currentSubpixelIsCentered = (j == 0);
                        for (int iteration = 0; iteration < iterationsToShow; iteration++) {
                            currentValue = LogisticFunction(currentR, currentValue);
                            if (currentValue >= xLowerBound && currentValue < xUpperBound) {
                                // normalize currentValue's range from [xLowerBound, xUpperBound) to [0, 1)
                                double normalizedCurrentValue =
                                        (currentValue - xLowerBound) / (xUpperBound - xLowerBound);
                                if (currentSubpixelIsCentered) {
                                    // sonification
                                    if (sonificationApplyInverseFourierTransform) {
                                        // bucket ranges from [0, ifftMaxFrequency - ifftMinFrequency)
                                        // audioDataIndex ranges from [ifftMinFrequency, ifftMaxFrequency - 1]
                                        double bucket = normalizedCurrentValue * (ifftMaxFrequency - ifftMinFrequency);
                                        int audioDataIndex = (int) std::floor(bucket + ifftMinFrequency);
                                        currentAudioData.at(audioDataIndex) += audioDataIncrement;
                                    } else {
                                        // normalizedCurrentValue ranges from [0, 1)
                                        // currentAudioDataPoint ranges from [-1, 1)
                                        double currentAudioDataPoint = normalizedCurrentValue * 2 - 1;
                                        // effectively divide the sample rate by rawAudioSampleRateFactor so that the oscillations between two values are not outside the range of human hearing
                                        for (int k = iteration * rawAudioSampleRateFactor;
                                             k < (iteration + 1) * rawAudioSampleRateFactor; k++) {
                                            currentAudioData.at(k) = currentAudioDataPoint;
                                        }
                                    }
                                    // graph
                                    if (doStartGraphAtIteration0) {
                                        currentGraphData.at(iteration + iterationsToSteadyState) = normalizedCurrentValue;
                                    } else {
                                        currentGraphData.at(iteration) = normalizedCurrentValue;
                                    }
                                    // plot
                                    currentDataPoints.at(iteration) = currentValue;
                                }
                                // subpixelHeight ranges from [0, canvasHeight)
                                // pixelHeight ranges from [0, canvasHeight - 1]
                                double subpixelHeight = normalizedCurrentValue * canvasHeight;
                                int pixelHeight = std::min((int) std::round(subpixelHeight), canvasHeight - 1);
                                // plot
                                if (currentSubpixelIsCentered) {
                                    currentPlotData.at(pixelHeight)++;
                                    if (subpixelHeight > maxY.at(i)) {
                                        maxY.at(i) = subpixelHeight;
                                    }
                                    if (subpixelHeight < minY.at(i)) {
                                        minY.at(i) = subpixelHeight;
                                    }
                                }
                                // visualization
                                if (antiAliasingEnabled) {
                                    // lowerPixelHeight and upperPixelHeight ranges from [0, canvasHeight - 1]
                                    int lowerPixelHeight = (int) std::floor(subpixelHeight);
                                    int upperPixelHeight = std::min((int) std::ceil(subpixelHeight), canvasHeight - 1);
                                    // won't check for lowerPixelHeight == upperPixelHeight since that will be very rare
                                    // subpixelHeight - lowerPixelHeight and upperPixelHeight - subpixelHeight might seem switched but they're not
                                    // TODO: remove lots of repetition by defining a function
                                    if (j < 0) {
                                        previousFrequencies.at(lowerPixelHeight) +=
                                                (upperPixelHeight - subpixelHeight) *
                                                (-j / (double) widthSamplesPerPixel) *
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
                                                (subpixelHeight - lowerPixelHeight) *
                                                (-j / (double) widthSamplesPerPixel) *
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
                                                (upperPixelHeight - subpixelHeight) *
                                                (j / (double) widthSamplesPerPixel) *
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
                                                (subpixelHeight - lowerPixelHeight) *
                                                (j / (double) widthSamplesPerPixel) *
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
                                    currentFrequencies.at(pixelHeight)++;
                                    if (currentMaxFrequency < currentFrequencies.at(pixelHeight)) {
                                        currentMaxFrequency = currentFrequencies.at(pixelHeight);
                                    }
                                }
                            }
                        }
                        iterationsCalculated++;
                    }
                }
                safeToRecalculate = true;
                parameterConditionVariable.notify_all();
                if (needsToRecalculate) {
                    std::cout << "cancelled calculation of logistic map\n";
                    needsToRecalculate = false;
                    return false;
                }
            }
        }
        iterationsCalculated = 0;
        numIterations = canvasWidth;
        calculationStage++;
        {
            // calculationStage 2: convert the calculations into formats ready for sonification and visualization
            int rangeRGBA[4];
            for (int i = 0; i < 4; i++) {
                rangeRGBA[i] = maxLogisticMapRGBA[i] - minLogisticMapRGBA[i];
            }
            double cutoffFrequency = std::exp(-logScalingFactor);
            for (int i = 0; i < canvasWidth; i++) {
                if (needsToRecalculate) {
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
                        for (int k = 0; k < 4; k++) {
                            imageData.at(pixelIndex + k) = backgroundRGBA[k];
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
                if (sonificationApplyInverseFourierTransform) {
                    {
                        const auto& currentAudioData = audioData.at(i);

                        // fftwIO is in halfcomplex format and currentAudioData has just the real components
                        // we want the imaginary parts to be 0 (for now at least)

                        // currentAudioData.resize(2 * ifftAudioSamples, 0.0);
                        // fftwIO = currentAudioData;
                        // NOTE: assumes size of fftwIn is 2 * ifftAudioSamples
                        for (int currentAudioIndex = 0; currentAudioIndex < 2 * ifftAudioSamples; currentAudioIndex++) {
                            if (currentAudioIndex < currentAudioData.size()) {
                                double magnitude = currentAudioData.at(currentAudioIndex);
                                double randomPhase = randomIfftPhase(rng);
                                fftwIn[currentAudioIndex][0] = magnitude * std::cos(randomPhase);
                                fftwIn[currentAudioIndex][1] = magnitude * std::sin(randomPhase);
                            } else {
                                fftwIn[currentAudioIndex][0] = 0.0;
                                fftwIn[currentAudioIndex][1] = 0.0;
                            }
                        }
                        fftw_execute(fftwPlan);
                        // audioData.at(i) = fftwIO;
                        audioData.at(i) = fftwOut;
                    }

                    // apply audio peak normalization
                    if (doIfftAudioPeakNormalization) {
                        auto& currentAudioData = audioData.at(i);
                        double maxAbsoluteAudio = 0;
                        // audioDataPoint is in [-1, 1)
                        for (int j = 0; j < currentAudioData.size(); j++) {
                            auto audioDataPoint = currentAudioData.at(j);
                            maxAbsoluteAudio = std::max(maxAbsoluteAudio, std::abs(audioDataPoint));
                        }
                        // for (const auto audioDataPoint : currentAudioData) {
                        //     maxAbsoluteAudio = std::max(maxAbsoluteAudio, std::abs(audioDataPoint - 0.5));
                        // }
                        // maxAbsoluteAudio is in [0, 0.5)
                        if (maxAbsoluteAudio != 1) {
                            for (auto& audioDataPoint : currentAudioData) {
                                audioDataPoint /= maxAbsoluteAudio;
                            }
                        }
                    }
                }

                iterationsCalculated++;
            }
            if (needsToRecalculate) {
                std::cout << "cancelled calculation of logistic map\n";
                needsToRecalculate = false;
                return false;
            }
        }

        std::cout << "calculated logistic map\n";
        currentlyCalculating = false;
        return true;
    }

public:

    void setStartingValue(double x0) {
        startingValue = x0;
    }
    
    void setDefaultEquationType(LogisticMap::DefaultEquationType equationType) {
        using enum LogisticMap::DefaultEquationType;
        switch (equationType) {
            case LOGISTIC:
                functionID = 0;

                rLowerBound = 2.75;
                rUpperBound = 4;
                xLowerBound = 0;
                xUpperBound = 1;
                break;
            case EXPONENTIAL:
                functionID = 1;

                rLowerBound = 1;
                rUpperBound = 30;
                xLowerBound = 0;
                xUpperBound = 11;
                break;
            case GAMMA:
                functionID = 2;

                rLowerBound = 1;
                rUpperBound = 9;
                xLowerBound = 0;
                xUpperBound = 10;
                break;
        }
    }
    
    void calculateLogisticMap() {
        bool calculationWasSuccessful;
        do {
            calculationWasSuccessful = tryCalculateLogisticMap();
            std::cout << calculationWasSuccessful << " success\n";
        }
        while (!calculationWasSuccessful);
    }
    void drawLogisticMap(emscripten::val canvas) {
        auto ctx = canvas.call<emscripten::val>("getContext", emscripten::val("2d"));
        auto canvasWidth = canvas["clientWidth"].as<int>();
        auto canvasHeight = canvas["clientHeight"].as<int>();

        std::cout << "drawing logistic map\n";
        // TODO: maybe make emscripten directly interpret a std::vector<char> as a Uint8ClampedArray
        ctx.call<void>("clearRect", emscripten::val(0), emscripten::val(0), canvas["width"], canvas["height"]);
        auto imageDataJsArray = emscripten::val(imageData);
        auto Uint8ClampedArray = emscripten::val::global("Uint8ClampedArray");
        auto imageDataJsUint8ClampedArray = Uint8ClampedArray.new_(imageDataJsArray);
        auto ImageData = emscripten::val::global("ImageData");
        auto imageDataJsObject = ImageData.new_(imageDataJsUint8ClampedArray, canvas["clientWidth"], canvas["clientHeight"]);
        ctx.call<void>("putImageData", imageDataJsObject, emscripten::val(0), emscripten::val(0));
        // no need to clear the canvas, this does it automatically
        std::cout << "drew logistic map\n";
    }

    void drawLogisticMapOverlay(emscripten::val canvas) {
        static std::vector<double> solidLinePatternArray = {};
        static emscripten::val solidLinePattern = emscripten::val(solidLinePatternArray);
        static std::vector<double> dashedLinePatternArray = {5, 5};
        static emscripten::val dashedLinePattern = emscripten::val(dashedLinePatternArray);

        emscripten::val ctx = canvas.call<emscripten::val>("getContext", emscripten::val("2d"));
        double currentX = currentMouseCoordinates[0];
        double currentY = currentMouseCoordinates[1];
        ctx.call<void>("clearRect", emscripten::val(0), emscripten::val(0), canvas["width"], canvas["height"]);
        if (currentlyCalculating) {
            return;
        }
        if (isCurrentlyPlaying) {
            // TODO: save and restore for all
            ctx.call<void>("save");
            ctx.set("strokeStyle", emscripten::val("rgba(0, 0, 0, 1)"));
            ctx.set("fillStyle", emscripten::val("rgba(0, 0, 0, 0.5)"));
            ctx.call<void>("fillRect", 0, canvasHeight - minY.at(currentXCoord), currentXCoord, -(maxY.at(currentXCoord) - minY.at(currentXCoord)));
            ctx.call<void>("strokeRect", 0, canvasHeight - minY.at(currentXCoord), currentXCoord, -(maxY.at(currentXCoord) - minY.at(currentXCoord)));
            ctx.call<void>("restore");
        } else {
            ctx.call<void>("setLineDash", dashedLinePattern);
            ctx.call<void>("beginPath");
            ctx.call<void>("moveTo", emscripten::val(currentXCoord), emscripten::val(0));
            ctx.call<void>("lineTo", emscripten::val(currentXCoord), canvas["height"]);
            ctx.call<void>("stroke");
        }
        if (enableOverlay and doShowOverlay) {
            ctx.call<void>("beginPath");
            ctx.call<void>("moveTo", emscripten::val(currentX), emscripten::val(0));
            ctx.call<void>("lineTo", emscripten::val(currentX), canvas["height"]);
            ctx.call<void>("stroke");
            ctx.call<void>("beginPath");
            ctx.call<void>("moveTo", emscripten::val(0), emscripten::val(currentY));
            ctx.call<void>("lineTo", canvas["width"], emscripten::val(currentY));
            ctx.call<void>("stroke");
            // NOTE: the label relies on the assumption that the canvases are centered
            std::string xText = "r = " + std::to_string(rLowerBound + (currentX / canvasWidth) * (rUpperBound - rLowerBound));
            std::string yText = "x = " + std::to_string(xLowerBound + (1 - currentY / canvasHeight) * (xUpperBound - xLowerBound));
            if (sonificationApplyInverseFourierTransform && isCurrentlyPlaying) {
                // remember that currentY starts from 0 at the top, not at the bottom
                yText = yText + " (" + std::to_string(ifftMinFrequency + (1 - currentY / canvasHeight) * (ifftMaxFrequency - ifftMinFrequency)) + " Hz)";
            }
            emscripten::val xTextMetrics = ctx.call<emscripten::val>("measureText", emscripten::val(xText));
            emscripten::val yTextMetrics = ctx.call<emscripten::val>("measureText", emscripten::val(yText));
            double xLeft    = xTextMetrics["actualBoundingBoxLeft"].as<double>();
            double xRight   = xTextMetrics["actualBoundingBoxRight"].as<double>();
            double xUp      = xTextMetrics["actualBoundingBoxAscent"].as<double>();
            double xDown    = xTextMetrics["actualBoundingBoxDescent"].as<double>();
            double yLeft    = yTextMetrics["actualBoundingBoxLeft"].as<double>();
            double yRight   = yTextMetrics["actualBoundingBoxRight"].as<double>();
            double yUp      = yTextMetrics["actualBoundingBoxAscent"].as<double>();
            double yDown    = yTextMetrics["actualBoundingBoxDescent"].as<double>();
            double greaterLeft = std::max(xLeft, yLeft);
            double greaterRight = std::max(xRight, yRight);
            double boundingBoxX = greaterLeft + greaterRight + 2 * overlayTextMargin;
            double boundingBoxY = xUp + 1.15 * (xDown + yUp) + yDown + 2 * overlayTextMargin;
            bool showLabel = true;
            std::array<int, 2> textboxLocation; // (1, 1) to put the textbox down and right, (-1, 1) to put it up and right, etc.
            if (boundingBoxX < canvasWidth - currentX) {
                textboxLocation[0] = 1;
            } else if (boundingBoxX < currentX) {
                textboxLocation[0] = -1;
            } else {
                showLabel = false;
            }
            if (boundingBoxY < canvasHeight - currentY) {
                textboxLocation[1] = 1;
            } else if (boundingBoxY < currentY) {
                textboxLocation[1] = -1;
            } else {
                showLabel = false;
            }
            if (showLabel) {
                ctx.set("fillStyle", emscripten::val("rgba(255, 255, 255, 0.5)"));
                ctx.call<void>("fillRect", currentX, currentY, textboxLocation[0] * boundingBoxX,
                               textboxLocation[1] * boundingBoxY);
                ctx.set("fillStyle", emscripten::val("black"));
                // NOTE: this part especially will break if the canvases are not centered
                ctx.call<void>("fillText", xText, currentX + textboxLocation[0] * boundingBoxX / 2, currentY + textboxLocation[1] * boundingBoxY / 2 - 1.15 * xDown);
                ctx.call<void>("fillText", yText, currentX + textboxLocation[0] * boundingBoxX / 2, currentY + textboxLocation[1] * boundingBoxY / 2 + 1.15 * yUp);
            }
        }
    }

    void drawPlot(emscripten::val canvas) {
        if (plotData.size() == 0 || plotData[0].size() == 0) {
            return;
        }
        emscripten::val ctx = canvas.call<emscripten::val>("getContext", emscripten::val("2d"));
        ctx.call<void>("clearRect", emscripten::val(0), emscripten::val(0), canvas["width"], canvas["height"]);
        auto& currentPlotData = plotData[currentXCoord];
        ctx.call<void>("beginPath");
        // NOTE: this assumes the plot canvas and the logistic map canvas have the same dimensions
        ctx.call<void>("moveTo", emscripten::val(canvasWidth), emscripten::val(canvasHeight));
        for (int i = 0; i < canvasHeight; i++) {
            ctx.call<void>("lineTo", emscripten::val(canvasWidth - currentPlotData[i]), emscripten::val(canvasHeight - i));
        }
        ctx.call<void>("stroke");
    }

    void drawWaveform(emscripten::val canvas) {
      // TODO: make a canvas class to interpret emscripten::val canvas
      if (plotData.size() == 0 || plotData[0].size() == 0) {
        return;
      }
      emscripten::val ctx = canvas.call<emscripten::val>("getContext", emscripten::val("2d"));
      // TODO: canvasWidth for drawPlot?
      auto canvasWidth = canvas["clientWidth"].as<int>();
      auto canvasHeight = canvas["clientHeight"].as<int>();
      ctx.call<void>("clearRect", emscripten::val(0), emscripten::val(0), canvas["width"], canvas["height"]);
      auto& currentPlotData = audioData.at(currentXCoord);
      ctx.call<void>("beginPath");
      // NOTE: this assumes the plot canvas and the logistic map canvas have the same dimensions
      ctx.call<void>("moveTo", emscripten::val(0), emscripten::val(0.5 * canvasHeight));
      for (int i = 0; i < std::min(currentPlotData.size(), (std::size_t) canvasWidth); i++) {
        ctx.call<void>("lineTo", emscripten::val(i), emscripten::val(0.5 * canvasHeight - 0.5 * canvasHeight * currentPlotData.at(i)));
      }
      ctx.call<void>("stroke");
    }

    void drawGraph(emscripten::val canvas) {
        static const double EXPANSION_FACTOR = 5;
        if (plotData.size() == 0 || plotData[0].size() == 0) {
            return;
        }
        emscripten::val ctx = canvas.call<emscripten::val>("getContext", emscripten::val("2d"));
        auto canvasWidth = canvas["clientWidth"].as<int>();
        auto canvasHeight = canvas["clientHeight"].as<int>();
        ctx.call<void>("clearRect", emscripten::val(0), emscripten::val(0), canvas["width"], canvas["height"]);
        auto& currentPlotData = graphData.at(currentXCoord);
        ctx.call<void>("beginPath");
        // NOTE: this assumes the plot canvas and the logistic map canvas have the same dimensions
        ctx.call<void>("moveTo", emscripten::val(0), emscripten::val(0.5 * canvasHeight));
        for (int i = 0; i < std::min(currentPlotData.size(), (std::size_t) std::ceil(canvasWidth / EXPANSION_FACTOR)); i++) {
            ctx.call<void>("lineTo", emscripten::val(i * EXPANSION_FACTOR), emscripten::val(canvasHeight * currentPlotData.at(i)));
        }
        ctx.call<void>("stroke");
    }

    void createGainNode(int nodeID) {
        auto& audioCtx = audioCtxWrapper.value();
        auto& currentGainNodeWrapper = gainNodes.at(nodeID);
        int audioBufferLength;
        if (sonificationApplyInverseFourierTransform) {
            audioBufferLength = ifftAudioSamples;
        } else {
            audioBufferLength = iterationsToShow * rawAudioSampleRateFactor;
        }
        auto audioBuffer = audioCtx.call<emscripten::val>("createBuffer", emscripten::val(1),
                                                          emscripten::val(audioBufferLength),
                                                          audioCtx["sampleRate"]);
        auto audioDataJsArray = emscripten::val(audioData.at(nodeID));
        auto Float32Array = emscripten::val::global("Float32Array");
        auto audioDataJsFloat32Array = Float32Array.new_(audioDataJsArray);
        audioBuffer.call<void>("copyToChannel", audioDataJsFloat32Array, emscripten::val(0));
        auto bufferSourceNode = audioCtx.call<emscripten::val>("createBufferSource");
        bufferSourceNode.set("loop", emscripten::val(true));
        bufferSourceNode.set("buffer", audioBuffer);
        emscripten::val GainNode = emscripten::val::global("GainNode");
        currentGainNodeWrapper = GainNode.new_(audioCtx);
        auto& gainNode = currentGainNodeWrapper.value();
        // NOTE: this is directly done with set instead of with setValueAtTime or else the cancelAndHoldAtTime from lerpPlayGainNode will cancel it
        gainNode["gain"].set("value", emscripten::val(0));
        bufferSourceNode.call<void>("start", emscripten::val(0));
        bufferSourceNode.call<void>("connect", gainNode);
        gainNode.call<void>("connect", audioCtx["destination"]);
    }

    void disconnectAllGainNodes() {
        for (auto& gainNodeOptWrapper : gainNodes) {
            if (gainNodeOptWrapper.has_value()) {
                gainNodeOptWrapper.value().call<void>("disconnect");
                gainNodeOptWrapper.reset();
            }
        }
    }

    void lerpPlayGainNode(int nodeID, double currentTime) {
        auto& currentGainNode = gainNodes.at(nodeID).value();
        currentGainNode["gain"].call<void>("cancelAndHoldAtTime", emscripten::val(currentTime));
        currentGainNode["gain"].call<void>("setValueAtTime", currentGainNode["gain"]["value"], emscripten::val(currentTime));
        currentGainNode["gain"].call<void>("linearRampToValueAtTime", emscripten::val(1), emscripten::val(currentTime + changeChannelAudioLerpTime));
    }

    void lerpStopGainNode(int nodeID, double currentTime) {
        auto& currentGainNode = gainNodes.at(nodeID).value();
        currentGainNode["gain"].call<void>("cancelAndHoldAtTime", emscripten::val(currentTime));
        currentGainNode["gain"].call<void>("setValueAtTime", currentGainNode["gain"]["value"], emscripten::val(currentTime));
        currentGainNode["gain"].call<void>("linearRampToValueAtTime", emscripten::val(0), emscripten::val(currentTime + changeChannelAudioLerpTime));
    }

    void sonifyLogisticMap() {
        if (!audioCtxWrapper.has_value()) {
            emscripten::val AudioContext = emscripten::val::global("AudioContext");
            audioCtxWrapper = AudioContext.new_();
        }
        if (!isCurrentlyPlaying) {
            auto document = emscripten::val::global("document");
            emscripten::val canvas = document.call<emscripten::val>("getElementById",
                                                                    emscripten::val("canvas-logistic-map-overlay"));
            emscripten::val ctx = canvas.call<emscripten::val>("getContext", emscripten::val("2d"));
            ctx.set("strokeStyle", emscripten::val("red"));
        }
        int x = currentMouseCoordinates.at(0);
        if (x < 0) {x = 0;}
        if (x >= canvasWidth) {x = canvasWidth - 1;}
        auto& currentGainNodeWrapper = gainNodes.at(x);
        if (!currentGainNodeWrapper.has_value()) {
            // This is a relatively intensive task so currentTime should be measured after this step
            createGainNode(x);
        }
        auto& audioCtx = audioCtxWrapper.value();
        double currentTime = audioCtx["currentTime"].as<double>();
        if (isCurrentlyPlaying) {
            lerpStopGainNode(currentXCoord, currentTime);
        } else {
            isCurrentlyPlaying = true;
        }
        currentXCoord = x;
        lerpPlayGainNode(x, currentTime);
    }

    void desonifyLogisticMap() {
        auto document = emscripten::val::global("document");
        emscripten::val canvas = document.call<emscripten::val>("getElementById", emscripten::val("canvas-logistic-map-overlay"));
        emscripten::val ctx = canvas.call<emscripten::val>("getContext", emscripten::val("2d"));
        ctx.set("strokeStyle", emscripten::val("black"));
        if (isCurrentlyPlaying) {
            auto& audioCtx = audioCtxWrapper.value();
            double currentTime = audioCtx["currentTime"].as<double>();
            lerpStopGainNode(currentXCoord, currentTime);
            isCurrentlyPlaying = false;
        }
    }

    void showOverlay() {
        doShowOverlay = true;
    }

    void hideOverlay() {
        doShowOverlay = false;
    }
};

LogisticMap& GetLogisticMap() {
    static auto logisticMap = LogisticMap();
    return logisticMap;
}

void PrintThrownExceptionToCerr(std::string prependText, const std::exception& e) {
    std::cerr << prependText << "\nMore Info:\n" << e.what() << "\n";
}

void PrintThrownNonExceptionToCerr(std::string prependText = "") {
    std::cerr << prependText << "Unknown non-exception thrown like exception.\n";
}

void SetInputToSuccess(emscripten::val input) {
    try {
        input["style"].set("backgroundColor", emscripten::val(""));
    } catch (const std::exception& e) {
        PrintThrownExceptionToCerr("Unknown exception.", e);
    } catch (...) {
        PrintThrownNonExceptionToCerr();
    }
}

void SetInputToFailure(emscripten::val input) {
  try {
    input["style"].set("backgroundColor", emscripten::val("pink"));
  } catch (const std::exception& e) {
    PrintThrownExceptionToCerr("Unknown exception.", e);
  } catch (...) {
    PrintThrownNonExceptionToCerr();
  }
}

bool ManipulateLogisticMap(bool reparameterizeLogisticMap, bool resizeLogisticMap, bool calculateLogisticMap, bool renderLogisticMap) {
    auto& logisticMap = GetLogisticMap();
    auto document = emscripten::val::global("document");
    auto canvas = document.call<emscripten::val>("getElementById", emscripten::val("canvas-logistic-map"));
    auto progress = document.call<emscripten::val>("getElementById", emscripten::val("progress-logistic-map"));
    auto progressbar = document.call<emscripten::val>("getElementById", emscripten::val("progressbar-logistic-map"));
    auto progressLabel = document.call<emscripten::val>("getElementById", emscripten::val("label-progressbar-logistic-map"));

    std::array<bool, 5> finishValues = {false, false, false, false, false};

    if (reparameterizeLogisticMap) {
        auto dropdown = document.call<emscripten::val>("getElementById", emscripten::val("equation"));
        std::string dropdownValue = dropdown["value"].as<std::string>();
        if (dropdownValue == "logistic") {
            logisticMap.setDefaultEquationType(LogisticMap::DefaultEquationType::LOGISTIC);
            finishValues.at(0) = true;
            SetInputToSuccess(dropdown);
        } else if (dropdownValue == "exponential") {
            logisticMap.setDefaultEquationType(LogisticMap::DefaultEquationType::EXPONENTIAL);
            finishValues.at(0) = true;
            SetInputToSuccess(dropdown);
        } else if (dropdownValue == "gamma") {
            logisticMap.setDefaultEquationType(LogisticMap::DefaultEquationType::GAMMA);
            finishValues.at(0) = true;
            SetInputToSuccess(dropdown);
        } else {
            std::cerr << "dropdown value of " + dropdownValue + " was not found.\n";
            finishValues.at(0) = false;
            SetInputToFailure(dropdown);
        }

        auto x0Input = document.call<emscripten::val>("getElementById", emscripten::val("x"));
        std::string x0InputVal = x0Input["value"].as<std::string>();
        double x0;
        try {
            x0 = std::stod(x0InputVal);
            SetInputToSuccess(x0Input);
        } catch (const std::invalid_argument& e) {
            std::cerr << "x0 was not a valid number: " + x0InputVal + "\nMore Info:\n" + e.what() + "\n";
            SetInputToFailure(x0Input);
        } catch (const std::out_of_range& e) {
            std::cerr << "x0 was out of range: " + x0InputVal + "\nMore Info:\n" + e.what() + "\n";
            SetInputToFailure(x0Input);
        } catch (const std::exception& e) {
            std::cerr << "Unknown exception.\nMore Info:\n" << e.what() << "\n";
            SetInputToFailure(x0Input);
        } catch (...) {
            std::cerr << "Unknown non-exception non-string thrown as exception.\n";
            SetInputToFailure(x0Input);
        }
        logisticMap.setStartingValue(x0);
    }

    // we subtract 2 * logisticMap.extraSymmetricalSamplesPerPixel because of the subpixel edges that go off the
    // centers of the edge pixels of the map, which we don't calculate
    if (resizeLogisticMap) {
        auto canvasWidth = canvas["clientWidth"].as<int>();
        auto canvasHeight = canvas["clientHeight"].as<int>();
        logisticMap.resizeLogisticMap(canvasWidth, canvasHeight);
        finishValues.at(1) = true;
    }
    if (calculateLogisticMap) {
        if (logisticMap.currentlyCalculating) {
            logisticMap.needsToRecalculate = true;
        } else {
            logisticMap.currentlyCalculating = true;
            // technically redundant since calculateLogisticMap will also set currentlyCalculating to true, but starting
            // the thread is slower so we will manually set this to be true in time for the renderLogisticMap logic
            std::thread calculations(&LogisticMap::calculateLogisticMap, &logisticMap); // equivalent to logisticMap.calculateLogisticMap()
            calculations.detach();
            progressbar["style"].set("visibility", emscripten::val("visible"));
        }
        finishValues.at(2) = true;
    }
    if (renderLogisticMap) {
        if (logisticMap.currentlyCalculating) {
            progress.set("value", emscripten::val(logisticMap.iterationsCalculated));
            progress.set("max", emscripten::val(logisticMap.numIterations));
            std::string calculationStageString = std::to_string(logisticMap.calculationStage) + "/" + std::to_string(logisticMap.maxCalculationStage);
            std::string iterationsCalculatedString = std::to_string(logisticMap.iterationsCalculated) + "/" + std::to_string(logisticMap.numIterations);
            std::string progressString = "Calculating Logistic Map:\nStage " + calculationStageString + "\nrValue " + iterationsCalculatedString;
            // TODO: use proper text node instead
            progressLabel.set("innerHTML", emscripten::val(progressString));
        } else {
            // disconnect gain nodes from previous calculations
            logisticMap.disconnectAllGainNodes();
            progressbar["style"].set("visibility", emscripten::val("hidden"));
            logisticMap.drawLogisticMap(canvas);
            finishValues.at(3) = true;
        }
    }
    return finishValues.at(0) && reparameterizeLogisticMap || finishValues.at(1) && resizeLogisticMap || finishValues.at(2) && calculateLogisticMap || finishValues.at(3) && renderLogisticMap;
}

void RenderLogisticMapOverlay(double DOMHighResTimeStamp)
{
    emscripten::val document = emscripten::val::global("document");
    emscripten::val canvas = document.call<emscripten::val>("getElementById", emscripten::val("canvas-logistic-map-overlay"));
    auto& logisticMap = GetLogisticMap();
    logisticMap.drawLogisticMapOverlay(canvas);
}

void RenderLogisticMap(double DOMHighResTimeStamp) {
    // TODO: only render once per frame
    bool manipulationWasFinished = ManipulateLogisticMap(false, false, false, true);
    emscripten::val window = emscripten::val::global("window");
    if (!manipulationWasFinished) {
        window.call<void>("requestAnimationFrame", emscripten::val::module_property("RenderLogisticMap"));
    }
}

void RenderPlot(double DOMHighResTimeStamp)
{
    emscripten::val document = emscripten::val::global("document");
    emscripten::val canvas = document.call<emscripten::val>("getElementById", emscripten::val("canvas-plot"));
    auto& logisticMap = GetLogisticMap();
    logisticMap.drawPlot(canvas);
}

void RenderWaveform(double DOMHighResTimeStamp)
{
  emscripten::val document = emscripten::val::global("document");
  emscripten::val canvas = document.call<emscripten::val>("getElementById", emscripten::val("canvas-waveform"));
  auto& logisticMap = GetLogisticMap();
  logisticMap.drawWaveform(canvas);
}

void RenderGraph(double DOMHighResTimeStamp)
{
    emscripten::val document = emscripten::val::global("document");
    emscripten::val canvas = document.call<emscripten::val>("getElementById", emscripten::val("canvas-graph"));
    auto& logisticMap = GetLogisticMap();
    logisticMap.drawGraph(canvas);
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
                window.call<void>("requestAnimationFrame", emscripten::val::module_property("RenderPlot"));
                window.call<void>("requestAnimationFrame", emscripten::val::module_property("RenderWaveform"));
                window.call<void>("requestAnimationFrame", emscripten::val::module_property("RenderGraph"));
            }
        }
    } else if (eventName == "mouseenter") {
        auto &logisticMap = GetLogisticMap();
        logisticMap.showOverlay();
    } else if (eventName == "mouseout") {
        auto& logisticMap = GetLogisticMap();
        logisticMap.hideOverlay();
    }
    emscripten::val window = emscripten::val::global("window");
    window.call<void>("requestAnimationFrame", emscripten::val::module_property("RenderLogisticMapOverlay"));
}

void InitializeCanvases(emscripten::val event)
{
    emscripten::val document = emscripten::val::global("document");

    auto progressLabel = document.call<emscripten::val>("getElementById", emscripten::val("label-progressbar-logistic-map"));
    std::string progressString = "Initializing...";
    // TODO: use proper text node instead
    progressLabel.set("innerHTML", emscripten::val(progressString));

    document.call<emscripten::val>("querySelectorAll", emscripten::val(".canvas")).call<void>("forEach", emscripten::val::module_property("InitializeCanvas"));
    emscripten::val window = emscripten::val::global("window");
    ManipulateLogisticMap(true, true, true, false);
    window.call<void>("requestAnimationFrame", emscripten::val::module_property("RenderPlot"));
    window.call<void>("requestAnimationFrame", emscripten::val::module_property("RenderLogisticMap"));
    window.call<void>("requestAnimationFrame", emscripten::val::module_property("RenderWaveform"));
    window.call<void>("requestAnimationFrame", emscripten::val::module_property("RenderGraph"));
    // TODO: make canvas-logistic-map-overlay not hardcoded
    auto logisticMapCanvasOverlay = document.call<emscripten::val>("getElementById", emscripten::val("canvas-logistic-map-overlay"));
    document.call<void>("addEventListener", emscripten::val("mousedown"), emscripten::val::module_property("InteractWithLogisticMapCanvas"));
    document.call<void>("addEventListener", emscripten::val("mouseup"), emscripten::val::module_property("InteractWithLogisticMapCanvas"));
    logisticMapCanvasOverlay.call<void>("addEventListener", emscripten::val("mousemove"), emscripten::val::module_property("InteractWithLogisticMapCanvas"));
    logisticMapCanvasOverlay.call<void>("addEventListener", emscripten::val("mouseenter"), emscripten::val::module_property("InteractWithLogisticMapCanvas"));
    logisticMapCanvasOverlay.call<void>("addEventListener", emscripten::val("mouseout"), emscripten::val::module_property("InteractWithLogisticMapCanvas"));
}

void InitializeAllSettings()
{
}

int main()
{
    emscripten::val window = emscripten::val::global("window");
    emscripten::val document = emscripten::val::global("document");
    window.call<void>("addEventListener", emscripten::val("resize"), emscripten::val::module_property("InitializeCanvases"));
    emscripten::val recalculateButton = document.call<emscripten::val>("getElementById", emscripten::val("recalculate-button"));
    recalculateButton.call<void>("addEventListener", emscripten::val("mouseup"), emscripten::val::module_property("InitializeCanvases"));
    /*
    document.call<void>("addEventListener", emscripten::val("mousedown"), emscripten::val::module_property("InteractWithCanvas"));
    document.call<void>("addEventListener", emscripten::val("mouseup"), emscripten::val::module_property("InteractWithCanvas"));
    document.call<void>("addEventListener", emscripten::val("mousemove"), emscripten::val::module_property("InteractWithCanvas"));
    */
    // initialize width and height of the canvas
    // TODO:
    //  window.call<void>("addEventListener", emscripten::val("onload"), emscripten::val::module_property("InitializeCanvases"));
    //  https://html.spec.whatwg.org/multipage/scripting.html#attr-script-async
    InitializeCanvases(emscripten::val::null());

    // retrieve all settings from localStorage and set the appropriate boxes to "checked" and put the appropriate data into preview
    InitializeAllSettings();

    return 0;
}

EMSCRIPTEN_BINDINGS(bindings)\
{\
  emscripten::function("InitializeCanvases", InitializeCanvases);\
  emscripten::function("InitializeCanvas", InitializeCanvas);\
  emscripten::function("RenderLogisticMapOverlay", RenderLogisticMapOverlay);\
  emscripten::function("RenderLogisticMap", RenderLogisticMap);\
  emscripten::function("RenderPlot", RenderPlot);\
  emscripten::function("RenderWaveform", RenderWaveform);\
  emscripten::function("RenderGraph", RenderGraph);\
  emscripten::function("InteractWithLogisticMapCanvas", InteractWithLogisticMapCanvas);\
};