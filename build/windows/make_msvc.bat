doxygen "C:\CSMP\CSMP_DOxygen\Doxyfile"
tar -cvf docbook.tar -C "C:\CSMP\CSMP_Documentation" docbook
tar -cvf html.tar -C "C:\CSMP\CSMP_Documentation" html
tar -cvf latex.tar -C "C:\CSMP\CSMP_Documentation" latex
tar -cvf man.tar -C "C:\CSMP\CSMP_Documentation" man

if exist .\CSMP_number_types.h (
	bjam --toolset=msvc-14.0 --build-type=complete --layout=system optimization=speed address-model=64 link=static threading=multi 
) else (
	bjam --toolset=msvc-14.0 --build-type=complete --layout=system optimization=speed address-model=64 link=static threading=multi 
	..\..\bin\windows\csmptype.exe
	bjam --toolset=msvc-14.0 --build-type=complete --layout=system optimization=speed address-model=64 link=static threading=multi 
)
copy .\CSMP_number_types.h ..\..\source\includes\model\

bjam --toolset=msvc-14.0 --build-type=complete --layout=system optimization=speed address-model=64 link=static threading=multi 

copy ..\..\bin\windows\*.lib ..\..\bin\data\
copy ..\..\bin\windows\*.dll ..\..\bin\data\
copy ..\..\bin\windows\csmptest.exe ..\..\bin\data\
cd ..\..\bin\data
.\csmptest.exe
