#include "ColorPalette_Test.h"

namespace csmp {

void ColorPalette_Test::run()
{
    ColorPalette cp(100,100);

    //Test Default variables assigned by the constructor
    _test(cp.PaletteSize()==256);
    _test(cp.Saturation()==1.0F);

    //Test PaletteSize
    cp.PaletteSize(250);
    _test(cp.PaletteSize()==250);

    //Test Blend
    cp.Blend(1.5F);
    _test(cp.Blend()==1.5F);

    //Test Lightness
    cp.Lightness(200.0F);
    _test(cp.Lightness()==200.0F);

    //Test Saturation
    cp.Saturation(1.3F);
    _test(cp.Saturation()==1.3F);


    //Testing the == operator
    ColorPalette cp_new;
    cp_new=cp;

    _test(cp_new.PaletteSize()==250);
    _test(cp_new.Blend()==1.5F);
    _test(cp_new.Lightness()==200.0F);
    _test(cp_new.Saturation()==1.3F);

    _test(cp_new.PaletteSize()==cp.PaletteSize());
    _test(cp_new.Blend()==cp.Blend());
    _test(cp_new.Lightness()==cp.Lightness());
    _test(cp_new.Saturation()==cp.Saturation());

    //Test CopyConstructor
    ColorPalette cp_newest(cp);
    _test(cp_newest.PaletteSize()==250);
    _test(cp_newest.Blend()==1.5F);
    _test(cp_newest.Lightness()==200.0F);
    _test(cp_newest.Saturation()==1.3F);

    _test(cp_newest.PaletteSize()==cp.PaletteSize());
    _test(cp_newest.Blend()==cp.Blend());
    _test(cp_newest.Lightness()==cp.Lightness());
    _test(cp_newest.Saturation()==cp.Saturation());

    ColorPalette cpw2file(0,256);
    cpw2file.PaletteSize(1000);
    cpw2file.MakeRainbowPalette();

    _test(cpw2file.WriteRgbColorPaletteFile("colorPaletteRGBTest.txt")==true);
    _test(cpw2file.WriteRgbColorPaletteFile()==true);

    cpw2file.WriteRgbColorPaletteToStdout(); //Manual TEST
    cpw2file.WriteRgbColorPaletteToStdoutWithLineNumbers(); //Manual TEST

    ColorPalette cpfromfile;
    cpfromfile.ReadRgbColorPaletteFile("colorPaletteRGBTest.txt");
    //_equal( cpfromfile.PaletteSize(), 1000, 0 );

    cp.ScaleColorRangeTo(0,100);

    // TODO: test size: std::cout << "PaletteSize:" << cpfromfile.PaletteSize();
   }
}
