#include <vector>
#include "TString.h"
#include "TSystem.h"
namespace hlfv {
std::vector<TString> getFileList(const TString &pattern)
{
    std::vector<TString> files;

    // 1) If pattern is a directory: list all *.root in it
    void *dirHandle = gSystem->OpenDirectory(pattern);
    if (dirHandle)
    {
        TString dirname = pattern;
        if (!dirname.EndsWith("/"))
            dirname += "/";

        const char *entry;
        while ((entry = gSystem->GetDirEntry(dirHandle)))
        {
            TString fname = entry;
            if (fname == "." || fname == "..") continue;
            if (!fname.EndsWith(".root")) continue;

            files.push_back(dirname + fname);
        }
        gSystem->FreeDirectory(dirHandle);
        return files;
    }

    // 2) If pattern is an explicit existing .root file
    if (pattern.EndsWith(".root") && !gSystem->AccessPathName(pattern))
    {
        files.push_back(pattern);
        return files;
    }

    // 3) Otherwise treat as pattern, e.g. /path/prefix*.root
    TString dirname = gSystem->DirName(pattern);
    TString basename = gSystem->BaseName(pattern);

    void *patDir = gSystem->OpenDirectory(dirname);
    if (!patDir) return files;

    // prefix before first '*'
    TString prefix = basename;
    Ssiz_t starPos = prefix.Index("*");
    if (starPos != kNPOS)
        prefix = prefix(0, starPos);

    const char *entry;
    while ((entry = gSystem->GetDirEntry(patDir)))
    {
        TString fname = entry;
        if (fname == "." || fname == "..") continue;
        if (!fname.BeginsWith(prefix) || !fname.EndsWith(".root")) continue;

        files.push_back(dirname + "/" + fname);
    }
    gSystem->FreeDirectory(patDir);

    return files;
}
}