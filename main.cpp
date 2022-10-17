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
    int iterationsToShow = 5000;
    int extraSymmetricalSamplesPerPixel = 2;
    double contrast = 25;
    double rLowerBound = 2.75;
    double rUpperBound = 3.99;
    double startingValue = 0.5;
    double xLowerBound = 0;
    double xUpperBound = 1;
    bool antiAliasingEnabled = true;

    auto document = emscripten::val::global("document");
    auto canvas = document.call<emscripten::val>("getElementById", emscripten::val("canvas-logistic-map"));
    auto ctx = canvas.call<emscripten::val>("getContext", emscripten::val("2d"));
    auto canvasWidth = canvas["clientWidth"].as<int>();
    auto canvasHeight = canvas["clientHeight"].as<int>();
    std::cout << "logistic map\n";

    double rPixelStep = (rUpperBound - rLowerBound) / canvasWidth;
    int widthSamplesPerPixel = 2 * extraSymmetricalSamplesPerPixel + 1;
    std::vector<double> columnFrequencies(canvasHeight, 0);
    std::vector<std::vector<double>> frequencies(canvasWidth, columnFrequencies);
    for (long int i = 0; i < canvasWidth; i++) {
        // TODO: make this async
        double pixelR = rLowerBound + rPixelStep * i;
        std::vector<double> *currentFrequencies = &frequencies.at(i);
        std::vector<double> *previousFrequencies;
        if (i > 0) {
            previousFrequencies = &frequencies.at(i - 1);
        } else {
            previousFrequencies = currentFrequencies;
        }
        std::vector<double> *nextFrequencies;
        if (i < canvasWidth - 1) {
            nextFrequencies = &frequencies.at(i + 1);
        } else {
            nextFrequencies = currentFrequencies;
        }
        for (int j = -extraSymmetricalSamplesPerPixel; j < extraSymmetricalSamplesPerPixel; j++) {
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
                    int scalingFactor;
                    if (antiAliasingEnabled) {
                        int lowerPixelHeight = std::floor(subpixelHeight);
                        int upperPixelHeight = std::ceil(subpixelHeight);
                        // won't check for lowerPixelHeight == upperPixelHeight since that will be very rare
                        if (j < 0) {
                            previousFrequencies->at(lowerPixelHeight) +=
                                    (subpixelHeight - lowerPixelHeight) * (-j / widthSamplesPerPixel);
                            currentFrequencies->at(lowerPixelHeight) += (subpixelHeight - lowerPixelHeight) *
                                                                        ((extraSymmetricalSamplesPerPixel + j) /
                                                                         widthSamplesPerPixel);
                            previousFrequencies->at(upperPixelHeight) +=
                                    (upperPixelHeight - subpixelHeight) * (-j / widthSamplesPerPixel);
                            currentFrequencies->at(upperPixelHeight) += (upperPixelHeight - subpixelHeight) *
                                                                        ((extraSymmetricalSamplesPerPixel + j) /
                                                                         widthSamplesPerPixel);
                        } else if (j == 0) {
                            currentFrequencies->at(lowerPixelHeight) += (subpixelHeight - lowerPixelHeight);
                            currentFrequencies->at(upperPixelHeight) += (upperPixelHeight - subpixelHeight);
                        } else {
                            currentFrequencies->at(lowerPixelHeight) +=
                                    (subpixelHeight - lowerPixelHeight) * (-j / widthSamplesPerPixel);
                            nextFrequencies->at(lowerPixelHeight) += (subpixelHeight - lowerPixelHeight) *
                                                                     ((extraSymmetricalSamplesPerPixel + j) /
                                                                      widthSamplesPerPixel);
                            currentFrequencies->at(upperPixelHeight) +=
                                    (upperPixelHeight - subpixelHeight) * (-j / widthSamplesPerPixel);
                            nextFrequencies->at(upperPixelHeight) += (upperPixelHeight - subpixelHeight) *
                                                                     ((extraSymmetricalSamplesPerPixel + j) /
                                                                      widthSamplesPerPixel);
                        }
                    } else {
                        int pixelHeight = std::round(subpixelHeight);
                        currentFrequencies->at(pixelHeight)++;
                    }
                }
            }
        }
    }

    std::vector<unsigned char> data(canvasWidth * canvasHeight * 4, 0);
    for (long int i = 0; i < canvasWidth; i++) {
        int pixelAlphaIndex = (canvasWidth * (canvasHeight - 1) + i) * 4 + 3;
        std::vector<double> currentFrequencies = frequencies.at(i);
        for (int j = 0; j < canvasHeight; j++) {
            double frequency = currentFrequencies.at(j);
            unsigned char shade = std::min(255, (int) std::round((255.0 * frequency)/(widthSamplesPerPixel * contrast)));
            if (shade != 0) {
                data.at(pixelAlphaIndex) = shade;
            }
            pixelAlphaIndex -= canvasWidth * 4;
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