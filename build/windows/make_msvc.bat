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
