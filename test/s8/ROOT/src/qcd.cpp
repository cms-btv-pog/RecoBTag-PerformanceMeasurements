#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>

#include <TCanvas.h>
#include <TClass.h>
#include <TFile.h>
#include <TH1.h>
#include <THStack.h>
#include <TKey.h>
#include <TList.h>
#include <TRint.h>

using std::cerr;
using std::cout;
using std::endl;
using std::string;

using boost::lexical_cast;

double luminosity = 0;

void MergeRootfile(TDirectory *target, TList *sourcelist );
void scale(TFile *file, TH1 *);

int main(int argc, char *argv[])
try
{
    // Test if sufficient number of arguments is specified.
    if (5 > argc)
        throw std::invalid_argument("usage: merge luminosity out.root in.root [more_in.root]");

    std::auto_ptr<TRint> application(new TRint("histInMemory", 0, 0));

    TFile *output = 0;
    TList *inputs = 0;
    try
    {
        ::luminosity = lexical_cast<double>(argv[1]);
        if (0 >= ::luminosity)
            throw std::runtime_error("Non-positive luminosity is specified");

        // Create output file
        output = new TFile(argv[2], "RECREATE");
        if (!output->IsOpen())
            throw std::runtime_error("Failed to open output file");

        // open specified input files to the List
        //
        inputs = new TList();
        cout << "Inputs" << endl;
        for( int i = 3; argc > i; ++i)
        {
            cout << " [+] " << argv[i] << endl;
            inputs->Add(TFile::Open(argv[i]));
        }

        // Call merge
        MergeRootfile(output, inputs);

        // memory cleanup
        delete inputs;
        delete output;
    }
    catch(const std::exception &error)
    {
        // memory cleanup in case of error
        if (inputs)
            delete inputs;

        if (output)
            delete output;

        throw;
    }

    return 0;
}
catch(const std::exception &error)
{
    cerr << error.what() << endl;

    return 1;
}

void MergeRootfile(TDirectory *target, TList *sourcelist )
{
    TString path( (char*)strstr( target->GetPath(), ":" ) );
    path.Remove( 0, 2 );

    TFile *first_source = (TFile*)sourcelist->First();
    first_source->cd( path );
    TDirectory *current_sourcedir = gDirectory;
    //gain time, do not add the objects in the list in memory
    Bool_t status = TH1::AddDirectoryStatus();
    TH1::AddDirectory(kFALSE);

    // loop over all keys in this directory
    TIter nextkey( current_sourcedir->GetListOfKeys() );
    TKey *key, *oldkey=0;
    while ( (key = (TKey*)nextkey()))
    {
        //keep only the highest cycle number for each key
        if (oldkey &&
            !strcmp(oldkey->GetName(),key->GetName()))

            continue;

        cout << "process key: " << key->GetName() << endl;

        // read object from first source file
        first_source->cd( path );
        TObject *obj = key->ReadObj();

        if (obj->IsA()->InheritsFrom(TH1::Class()))
        {
            // descendant of TH1 -> merge it
            //
            TH1 *h1 = (TH1*) obj;
            scale(first_source, h1);

            // loop over all source files and add the content of the
            // correspondant histogram to the one pointed to by "h1"
            for(TFile *nextsource = (TFile*)sourcelist->After(first_source);
                nextsource;
                nextsource = (TFile*)sourcelist->After(nextsource))
            {
                // make sure we are at the correct directory level by cd'ing to path
                //
                nextsource->cd( path );
                TKey *key2 = (TKey*)gDirectory->GetListOfKeys()->FindObject(h1->GetName());
                if (key2)
                {
                    TH1 *h2 = (TH1*)key2->ReadObj();
                    scale(nextsource, h2);

                    h1->Add(h2);

                    delete h2;
                }
            }
        }
        else if ( obj->IsA()->InheritsFrom( TDirectory::Class() ) )
        {
            // it's a subdirectory

            cout << "Found subdirectory " << obj->GetName() << endl;

            // create a new subdir of same name and title in the target file
            target->cd();
            TDirectory *newdir = target->mkdir( obj->GetName(), obj->GetTitle() );

            // newdir is now the starting point of another round of merging
            // newdir still knows its depth within the target file via
            // GetPath(), so we can still figure out where we are in the recursion
            MergeRootfile( newdir, sourcelist );

        }
        else
        {

            // object is of no type that we know or can handle
            cout << "Unknown object type, name: "
            << obj->GetName() << " title: " << obj->GetTitle() << endl;

            continue;
        }

        // now write the merged histogram (which is "in" obj) to the target file
        // note that this will just store obj in the current directory level,
        // which is not persistent until the complete directory itself is stored
        // by "target->Write()" below
        if ( obj ) {
            target->cd();

            //!!if the object is a tree, it is stored in globChain...
            obj->Write( key->GetName() );
        }

    } // while ( ( TKey *key = (TKey*)nextkey() ) )

    // save modifications to target file
    target->SaveSelf(kTRUE);
    TH1::AddDirectory(status);
}

void scale(TFile *file, TH1 *hist)
{
    boost::filesystem::path path(file->GetName());
    std::string filename = boost::filesystem::basename(path);

    if ("pt15to20" == filename)
    {
        hist->Scale(579200000 * 0.00254 / 2884915 * luminosity);
    }
    else if ("pt20to30" == filename)
    {
        hist->Scale(236300000 * 0.00518 / 11461085 * luminosity);
    }
    else if ("pt30to50" == filename)
    {
        hist->Scale(53070000 * 0.01090 / 11431864 * luminosity);
    }
    else if ("pt50to80" == filename)
    {
        hist->Scale(6351000 * 0.02274 / 10748755 * luminosity);
    }
    else if ("pt80to120" == filename)
    {
        hist->Scale(785100 * 0.037 / 3191979 * luminosity);
    }
    else if ("pt120to150" == filename)
    {
        hist->Scale(92950 * 0.04777 / 998503 * luminosity);
    }
    else if ("pt150" == filename)
    {
        hist->Scale(47580 * 0.05964 / 1022541 * luminosity);
    }
    
    else
        cerr << "Didn't understand filename: " << filename
            << " Input is skipped" << endl;
}
