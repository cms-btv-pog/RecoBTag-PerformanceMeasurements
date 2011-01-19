#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <TCanvas.h>
#include <TClass.h>
#include <TFile.h>
#include <TH1.h>
#include <THStack.h>
#include <TKey.h>
#include <TList.h>
#include <TRint.h>
#include <TLegend.h>

using std::cerr;
using std::cout;
using std::endl;
using std::string;

double luminosity = 31.6;

class HeapObject
{
    public:
        HeapObject()
        {
            stack = 0;
            canvas = 0;
            legend = 0;
        }

        THStack *stack;
        TCanvas *canvas;
        TLegend *legend;
};

std::vector<HeapObject> heaps;

void MergeRootfile(TDirectory *target, TList *sourcelist );
std::string style(TFile *file, TH1 *);
TLegend *createLegend(const string &name = "");

int main(int argc, char *argv[])
try
{
    // Test if sufficient number of arguments is specified.
    if (4 > argc)
        throw std::invalid_argument("usage: merge out.root in.root [more_in.root]");

    std::auto_ptr<TRint> application(new TRint("histInMemory", 0, 0));

    TFile *output = 0;
    TList *inputs = 0;
    try
    {
        // Create output file
        output = new TFile(argv[1], "RECREATE");
        if (!output->IsOpen())
            throw std::runtime_error("Failed to open output file");

        // open specified input files to the List
        //
        inputs = new TList();
        cout << "Inputs" << endl;
        for( int i = 2; argc > i; ++i)
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

    application->Run(true);

    for(std::vector<HeapObject>::iterator heap = heaps.begin();
        heaps.end() != heap;
        ++heap)
    {
        if (heap->stack)
        {
            TIter next(heap->stack->GetHists());
            while(TObject *obj = next())
                delete obj;

            delete heap->stack;
        }

        if (heap->canvas)
            delete heap->canvas;

        if (heap->legend)
            delete heap->legend;
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
            THStack *stack = new THStack();
            HeapObject heap;
            heap.stack = stack;

            TLegend *legend = createLegend();
            heap.legend = legend;

            TH1 *h1 = (TH1*) obj;
            std::string label = style(first_source, h1);
            if (!label.empty())
            {
                legend->AddEntry(h1, label.c_str());
                stack->Add(h1);
            }

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
                    label = style(nextsource, h2);
                    if (!label.empty())
                    {
                        TH1 *clone = (TH1*) h2->Clone();
                        legend->AddEntry(clone, label.c_str());
                        stack->Add(clone);
                    }

                    delete h2;
                }
            }

            TCanvas *canvas = new TCanvas();
            heap.canvas = canvas;

            heaps.push_back(heap);
            
            //canvas->SetWindowSize(640, 480);
            stack->Draw("nostack e1");
            legend->Draw();

            canvas->RedrawAxis();
            canvas->Update();

            obj = canvas;
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

std::string style(TFile *file, TH1 *hist)
{
    std::string filename = file->GetName();

    if ("data.root" == filename)
    {
        hist->SetLineWidth(2);

        return "data";
    }

    if ("ppmux.root" == filename)
    {
        hist->SetLineWidth(4);
        hist->SetLineColor(kAzure + 1);

        hist->Scale(48440000 * 1.76 / 886239 * luminosity);

        return "ppmux";
    }

    if ("ttbar.root" == filename)
    {
        hist->SetLineWidth(4);
        hist->SetLineColor(kRed);

        hist->Scale(157.5 / 996859 * luminosity);

        return "ttbar";
    }

    if ("pt15.root" == filename)
    {
        hist->SetLineColor(kGray+1);
        hist->SetLineWidth(2);

        hist->Scale(874100000 * 0.0039 / 2986966 * luminosity);

        return "pt15";
    }

    if ("pt30.root" == filename)
    {
        hist->SetLineColor(kGreen + 1);
        hist->SetLineWidth(2);

        hist->Scale(61160000 * 0.0126 / 7628669 * luminosity);

        return "pt30";
    }

    if ("pt50.root" == filename)
    {
        hist->SetLineColor(kViolet);
        hist->SetLineWidth(2);

        hist->Scale(7290000 * 0.0246 / 4262850 * luminosity);

        return "pt150";
    }

    if ("pt150.root" == filename)
    {
        hist->SetLineColor(kBlue);
        hist->SetLineWidth(2);

        hist->Scale(48100 * 0.056 / 489271 * luminosity);

        return "pt150";
    }

    cerr << "Didn't understand filename: " << filename
        << " Input is not used" << endl;

    return "";
}

TLegend *createLegend(const string &name)
{
    TLegend *legend = new TLegend( .6, .6, .8, .9);
    if (!name.empty())
        legend->SetHeader(name.c_str());

    legend->SetMargin(0.12);
    legend->SetTextSize(0.035);
    legend->SetFillColor(10);
    legend->SetBorderSize(0);

    return legend;
}

