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


double LogisticFunction(double r, double input) {
    return r * input * (1 - input);
}

void RenderLogisticMap(double DOMHighResTimeStamp)
{
    int iterationsToSteadyState = 1000;
    int iterationsToShow = 2000;
    int extraSymmetricalSamplesPerPixel = 2;
    double logFactor = 8;
    double rLowerBound = 2.75;
    double rUpperBound = 3.99;
    double startingValue = 0.5;
    double xLowerBound = 0;
    double xUpperBound = 1;
    bool antiAliasingEnabled = true;
    bool logarithmicShadingEnabled = true; // alternative is linear shading
    unsigned char backgroundRGBA[4] = {255, 255, 255, 0};
    unsigned char maxLogisticMapRGBA[4] = {0, 0, 0, 255};
    unsigned char minLogisticMapRGBA[4] = {255, 255, 255, 0};

    auto document = emscripten::val::global("document");
    auto canvas = document.call<emscripten::val>("getElementById", emscripten::val("canvas-logistic-map"));
    auto ctx = canvas.call<emscripten::val>("getContext", emscripten::val("2d"));
    auto canvasWidth = canvas["clientWidth"].as<int>();
    auto canvasHeight = canvas["clientHeight"].as<int>();

    std::cout << "calculating logistic map\n";

    double rPixelStep = (rUpperBound - rLowerBound) / canvasWidth;
    int widthSamplesPerPixel = 2 * extraSymmetricalSamplesPerPixel + 1;
    std::vector<double> columnFrequencies(canvasHeight, 0);
    std::vector<std::vector<double>> frequencies(canvasWidth, columnFrequencies);
    std::vector<double> maxFrequencies(canvasWidth, 0);
    double scalingFactor = 1.0 / (widthSamplesPerPixel * iterationsToShow);
    for (long int i = 0; i < canvasWidth; i++) {
        // TODO: make this async
        double pixelR = rLowerBound + rPixelStep * i;
        auto& currentFrequencies = frequencies.at(i);
        auto& currentMaxFrequency = maxFrequencies.at(i);
        int previousFrequenciesIndex;
        if (i != 0) {
            previousFrequenciesIndex = i - 1;
        } else {
            previousFrequenciesIndex = i;
        }
        auto& previousFrequencies = frequencies.at(previousFrequenciesIndex);
        auto& previousMaxFrequency = maxFrequencies.at(previousFrequenciesIndex);
        int nextFrequenciesIndex;
        if (i != canvasWidth - 1) {
            nextFrequenciesIndex = i + 1;
        } else {
            nextFrequenciesIndex = i;
        }
        auto& nextFrequencies = frequencies.at(nextFrequenciesIndex);
        auto& nextMaxFrequency = maxFrequencies.at(nextFrequenciesIndex);
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
                    double subpixelHeight = canvasHeight * (currentValue - xLowerBound) / (xUpperBound - xLowerBound);
                    if (antiAliasingEnabled) {
                        int lowerPixelHeight = std::floor(subpixelHeight);
                        int upperPixelHeight = std::ceil(subpixelHeight);
                        // won't check for lowerPixelHeight == upperPixelHeight since that will be very rare
                        // subpixelHeight - lowerPixelHeight and upperPixelHeight - subpixelHeight might seem switched but they're not
                        // TODO: remove lots of repetition by defining a function
                        if (j < 0) {
                            previousFrequencies.at(lowerPixelHeight) +=
                                    (upperPixelHeight - subpixelHeight) * (-j / (double) widthSamplesPerPixel) * scalingFactor;
                            if (previousMaxFrequency < previousFrequencies.at(lowerPixelHeight)) {
                                previousMaxFrequency = previousFrequencies.at(lowerPixelHeight);
                            }
                            currentFrequencies.at(lowerPixelHeight) += (upperPixelHeight - subpixelHeight) *
                                    ((widthSamplesPerPixel + j) / (double) widthSamplesPerPixel) * scalingFactor;
                            if (currentMaxFrequency < currentFrequencies.at(lowerPixelHeight)) {
                                currentMaxFrequency = currentFrequencies.at(lowerPixelHeight);
                            }
                            previousFrequencies.at(upperPixelHeight) +=
                                    (subpixelHeight - lowerPixelHeight) * (-j / (double) widthSamplesPerPixel) * scalingFactor;
                            if (previousMaxFrequency < previousFrequencies.at(upperPixelHeight)) {
                                previousMaxFrequency = previousFrequencies.at(upperPixelHeight);
                            }
                            currentFrequencies.at(upperPixelHeight) += (subpixelHeight - lowerPixelHeight) *
                                                                        ((widthSamplesPerPixel + j) / (double) widthSamplesPerPixel) * scalingFactor;
                            if (currentMaxFrequency < currentFrequencies.at(upperPixelHeight)) {
                                currentMaxFrequency = currentFrequencies.at(upperPixelHeight);
                            }
                        } else if (j == 0) {
                            currentFrequencies.at(lowerPixelHeight) += (upperPixelHeight - subpixelHeight) * scalingFactor;
                            if (currentMaxFrequency < currentFrequencies.at(lowerPixelHeight)) {
                                currentMaxFrequency = currentFrequencies.at(lowerPixelHeight);
                            }
                            currentFrequencies.at(upperPixelHeight) += (subpixelHeight - lowerPixelHeight) * scalingFactor;
                            if (currentMaxFrequency < currentFrequencies.at(upperPixelHeight)) {
                                currentMaxFrequency = currentFrequencies.at(upperPixelHeight);
                            }
                        } else {
                            currentFrequencies.at(lowerPixelHeight) +=
                                    (upperPixelHeight - subpixelHeight) * (j / (double) widthSamplesPerPixel) * scalingFactor;
                            if (currentMaxFrequency < currentFrequencies.at(lowerPixelHeight)) {
                                currentMaxFrequency = currentFrequencies.at(lowerPixelHeight);
                            }
                            nextFrequencies.at(lowerPixelHeight) += (upperPixelHeight - subpixelHeight) *
                                    ((widthSamplesPerPixel - j) / (double) widthSamplesPerPixel) * scalingFactor;
                            if (nextMaxFrequency < nextFrequencies.at(lowerPixelHeight)) {
                                nextMaxFrequency = nextFrequencies.at(lowerPixelHeight);
                            }
                            currentFrequencies.at(upperPixelHeight) +=
                                    (subpixelHeight - lowerPixelHeight) * (j / (double) widthSamplesPerPixel) * scalingFactor;
                            if (currentMaxFrequency < currentFrequencies.at(upperPixelHeight)) {
                                currentMaxFrequency = currentFrequencies.at(upperPixelHeight);
                            }
                            nextFrequencies.at(upperPixelHeight) += (subpixelHeight - lowerPixelHeight) *
                                    ((widthSamplesPerPixel - j) / (double) widthSamplesPerPixel) * scalingFactor;
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
        }
    }

    std::cout << "drawing logistic map\n";

    int rangeRGBA[4];
    for (int i = 0; i < 4; i++) {
        rangeRGBA[i] = maxLogisticMapRGBA[i] - minLogisticMapRGBA[i];
    }
    std::vector<unsigned char> data(canvasWidth * canvasHeight * 4, 0);
    double generalScalingFactor;
    if (logarithmicShadingEnabled) {
        generalScalingFactor = -std::log(widthSamplesPerPixel) / logFactor;
        // generalScalingFactor = std::log(1 / widthSamplesPerPixel) / logFactor
    } else {
        generalScalingFactor = 1.0 / widthSamplesPerPixel;
    }
    for (long int i = 0; i < canvasWidth; i++) {
        int pixelIndex = (canvasWidth * (canvasHeight - 1) + i) * 4;
        std::vector<double> currentFrequencies = frequencies.at(i);
        double max = maxFrequencies.at(i);
        double columnScalingFactor;
        if (logarithmicShadingEnabled) {
            columnScalingFactor = generalScalingFactor - std::log(max) / logFactor;
            // columnScalingFactor = std::log(1 / (widthSamplesPerPixel * max)) / logFactor
        } else {
            columnScalingFactor = generalScalingFactor / max;
            // columnScalingFactor = 1 / (max * widthSamplesPerPixel)
        }
        for (int j = 0; j < canvasHeight; j++) {
            auto currentFrequency = currentFrequencies.at(j);
            if (currentFrequency == 0) {
                for (int i = 0; i < 4; i++) {
                    data.at(pixelIndex + i) = backgroundRGBA[i];
                }
            } else {
                double clampedPixelScalingFactor;
                if (logarithmicShadingEnabled) {
                    double pixelScalingFactor = std::log(currentFrequency) / logFactor + columnScalingFactor;
                    // pixelScalingFactor = std::log(currentFrequency / (widthSamplesPerPixel * max)) / logFactor
                    clampedPixelScalingFactor = std::min(255.0, std::max(0.0, pixelScalingFactor + 1));
                    // pixelScalingFactor scales from [-1, 0] for x in [e^-logFactor, 1] and [-infinity, 0] for x in [0, 1] barring anti-aliasing inaccuracy
                    // clampedPixelScalingFactor scales from [0, 1] for x in [e^-logFactor, 1] and for x in [0, 1]
                } else {
                    double pixelScalingFactor = currentFrequency * columnScalingFactor;
                    clampedPixelScalingFactor = std::min(255.0, std::max(0.0, pixelScalingFactor));
                    // pixelScalingFactor = currentFrequency / (max * widthSamplesPerPixel)
                    // pixelScalingFactor scales from [0, 1] for x in [0, 1] barring anti-aliasing inaccuracy
                }
                for (int i = 0; i < 4; i++) {
                    data.at(pixelIndex + i) = clampedPixelScalingFactor * rangeRGBA[i] + minLogisticMapRGBA[i];
                }
            }
            /*
            double pixelShade; // 0 to 1
            if (logarithmicShadingEnabled) {
                double logarithmicFrequency = std::log(logFactor * (currentFrequencies.at(j) / widthSamplesPerPixel) + 1) / (std::log(logFactor + 1) * max);
                shade = std::max(0, std::min(255, (int) std::round(255 * logarithmicFrequency)));
            } else {
                shade = std::min(255, (int) std::round(255 * currentFrequencies.at(j) / (widthSamplesPerPixel * max)));
            }
            if (shade != 0) {
                // this is so that background pixels are different from shaded pixels in programs that don't accept alpha values
                data.at(pixelIndex) = logisticMapRGBA[0];
                data.at(pixelIndex + 1) = 0;
                data.at(pixelIndex + 2) = 0;
                data.at(pixelIndex + 3) = shade;
            }
            */
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
    std::cout << "logistic map finished\n";
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