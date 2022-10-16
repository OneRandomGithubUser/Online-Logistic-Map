#include <emscripten/val.h>
#include <emscripten/bind.h>
#include <emscripten/emscripten.h>
#include <array>
#include <string>
#include <vector>
#include <string>
#include <iostream>

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


void Render(emscripten::val canvas)
{
    emscripten::val ctx = canvas.call<emscripten::val>("getContext", emscripten::val("2d"));
}

void RenderLogisticMap(double DOMHighResTimeStamp)
{
    int iterationsToSteadyState = 100;
    int iterationsToShow = 100;
    int widthSamplesPerPixel = 1;
    double rLowerBound = 0;
    double rUpperBound = 4;

    auto document = emscripten::val::global("document");
    auto canvas = document.call<emscripten::val>("getElementById", emscripten::val("canvas-logistic-map"));
    auto ctx = canvas.call<emscripten::val>("getContext", emscripten::val("2d"));
    auto canvasWidth = canvas["clientWidth"].as<int>();
    auto canvasHeight = canvas["clientHeight"].as<int>();
    std::cout << "logistic map\n";
    std::vector<unsigned char> data; // TODO: maybe make emscripten directly interpret a std::vector<char> as a Uint8ClampedArray
    data.reserve(canvasWidth * canvasHeight * 4);
    for (long int i = 0; i < canvasWidth * canvasHeight; i++) {
        for (int j = 0; j < widthSamplesPerPixel; j++)
        {
            //
        }
        // imageDataVal["data"].set(emscripten::val(i), emscripten::val(127));
        data.emplace_back(128);
        data.emplace_back(128);
        data.emplace_back(128);
        data.emplace_back(128);
    }
    std::cout << "1\n";
    auto data2 = emscripten::val(data);
    std::cout << "2\n";
    auto Uint8ClampedArray = emscripten::val::global("Uint8ClampedArray");
    auto data4 = Uint8ClampedArray.new_(data2);
    std::cout << "4\n";
    auto ImageData = emscripten::val::global("ImageData");
    auto data6 = ImageData.new_(data4, canvas["clientWidth"], canvas["clientHeight"]);
    std::cout << "6\n";
    ctx.call<void>("putImageData", data6, emscripten::val(0), emscripten::val(0));
    std::cout << "7\n";

    // ctx.call<void>("putImageData", imageDataVal, emscripten::val(0), emscripten::val(0));
    //for (int i = 0; i < widthSamplesPerPixel * canvasWidth; i++) {}
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
    window.call<void>("requestAnimationFrame", emscripten::val::module_property("RenderPlot"));
    window.call<void>("requestAnimationFrame", emscripten::val::module_property("RenderLogisticMap"));
}

void InitializeCanvases(emscripten::val event)
{
    emscripten::val document = emscripten::val::global("document");
    document.call<emscripten::val>("querySelectorAll", emscripten::val(".canvas")).call<void>("forEach", emscripten::val::module_property("InitializeCanvas"));
}

void InteractWithCanvas(emscripten::val event)
{
    emscripten::val window = emscripten::val::global("window");
    window.call<void>("requestAnimationFrame", emscripten::val::module_property("RenderPlot"));
    window.call<void>("requestAnimationFrame", emscripten::val::module_property("RenderLogisticMap"));
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
    document.call<emscripten::val>("querySelectorAll", emscripten::val(".canvas")).call<void>("forEach", emscripten::val::module_property("InitializeCanvas"));

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
  emscripten::function("Render", Render);\
};