#include "../include/utils.h"

#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TList.h>
#include <TSystem.h>

namespace hlfv {

std::vector<TString> getFileList(const TString &pattern)
{
    std::vector<TString> files;
    TString dirname = gSystem->DirName(pattern);
    TString basename = gSystem->BaseName(pattern);

    TSystemDirectory dir(dirname, dirname);
    TList *list = dir.GetListOfFiles();
    if (!list) return files;

    TIter next(list);
    TSystemFile *file;
    while ((file = (TSystemFile *)next()))
    {
        TString fname = file->GetName();
        if (!file->IsDirectory() && fname.Contains(basename.ReplaceAll("*", "")))
            files.push_back(dirname + "/" + fname);
    }
    delete list;
    return files;
}

} // namespace hlfv
