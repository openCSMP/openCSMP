// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "OS_Utilities.h"
#include "Exception.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <time.h>
#include <iostream>
#ifdef _WIN32
#include <windows.h>
#include <iostream>
#include <sstream>
#include <iomanip>
#else
#include <unistd.h>
#endif

using namespace std;

namespace csmp {

std::string currentDirectorySymbol()
{
#ifdef _WIN32
    return std::string(".\\");
#else
    return std::string("./");
#endif
}

std::string directorySymbol()
{
#ifdef _WIN32
    return std::string("\\");
#else
    return std::string("/");
#endif
}

void createDirectoryIfDoesntExist(  const std::string& directory_name )
{
#ifdef _WIN32
    CreateDirectory(directory_name.c_str(), NULL);
#else
    struct stat st = {0};
    if (stat(directory_name.c_str(), &st) == -1)
    {
        mkdir(directory_name.c_str(), 0700);
    }
#endif
}

std::string CompareFileModifiedTimeStamps(std::string file0_name,std::string file1_name){
    std::string earlier;
#ifndef _WIN32
    std::string mod_date0=getFileModificationTime((file0_name).c_str());
    std::string mod_date1=getFileModificationTime((file1_name).c_str());
    std::cout<<" mod date: "<<file0_name<<" "<<mod_date0<<endl;
    std::cout<<" mod date: "<<file1_name<<" "<<mod_date1<<endl;
    time_t t0, t1;
    struct tm tm0, tm1;

    //(1) convert `String to tm`:
    if(strptime(mod_date0.c_str(), "%X-%x",&tm0) == NULL)
        std::cout<<"strptime failed"<<endl;
    if(strptime(mod_date1.c_str(), "%X-%x",&tm1) == NULL)
        std::cout<<"strptime failed"<<endl;

    //(2)   convert `tm to time_t`:
    t0 = mktime(&tm0);
    t1 = mktime(&tm1);

    double seconds = difftime(t0,t1);
    if (seconds> 0)
        earlier=file1_name;
    else
        earlier=file0_name;
#else
    HANDLE hFile1,hFile0;
    hFile0 = CreateFile(file0_name.c_str(), GENERIC_READ, FILE_SHARE_READ, NULL,
                       OPEN_EXISTING, 0, NULL);

    hFile1 = CreateFile(file1_name.c_str(), GENERIC_READ, FILE_SHARE_READ, NULL,
                       OPEN_EXISTING, 0, NULL);

    FILETIME creationTime0,creationTime1,
             lpLastAccessTime0,lpLastAccessTime1,
             lastWriteTime0,lastWriteTime1;
    bool err = GetFileTime( hFile0, &creationTime0, &lpLastAccessTime0, &lastWriteTime0 );
    err = GetFileTime( hFile1, &creationTime1, &lpLastAccessTime1, &lastWriteTime1 );

    if (CompareFileTime(&lastWriteTime0,&lastWriteTime1)>0)  // a -1 implies file0 was written earlier (file1 is newer),
        earlier=file0_name;                                  // while +1 means file0 was written later (file1 is older).
    else
        earlier=file1_name;
    cout<<"restart file comparison (creation)  : "<<CompareFileTime(&creationTime0,&creationTime1)<<endl;
    cout<<"restart file comparison (access)    : "<<CompareFileTime(&lpLastAccessTime0,&lpLastAccessTime1)<<endl;
    cout<<"restart file comparison (write time): "<<CompareFileTime(&lastWriteTime0,&lastWriteTime1)<<endl;
    //    cin.get();
    CloseHandle(hFile0);
    CloseHandle(hFile1);
#endif
    return earlier;
}

std::string getFileModificationTime(const char *filePath)
{
#ifdef _WIN32
    throw csmp::Exception(FATAL_ERROR,"getFileModificationTime","getting file modification time is not available in Windows. Please write it!");
#endif
    struct stat attrib;
    stat(filePath, &attrib);
    char date[10];
    //strftime(date, 20, "%d-%m-%y", localtime(&(attrib.st_ctime)));
    strftime(date, 20, "%X-%x", localtime(&(attrib.st_mtime)));
    std::string date_str(date);
    date[0] = '\0';
    return date_str;
}

}
 // end namespace csmp
