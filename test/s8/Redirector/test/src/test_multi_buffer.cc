/**
 * test_multi_buffer
 *
 *
 * Created by Samvel Khalatian on Sep 10, 2010
 * Copyright 2010, All rights reserved
 */

#include <iostream>
#include <fstream>

#include "interface/Debug.h"
#include "interface/MultiStreamBuffer.h"
#include "interface/Redirector.h"

using std::cerr;
using std::cout;
using std::clog;
using std::endl;
using std::streambuf;

using top::Debug;
using top::MultiStreamBuffer;
using top::Redirector;

void test1();
void test2();
void test3();

int main(int argc, char *argv[])
{
    cout << "Test1" << endl;
    cout << "-----" << endl;
    cout << "  Everything is output to: COUT and CERR." << endl;
    cout << endl;

    test1();

    cout << endl;
    cout << "Test2" << endl;
    cout << "-----" << endl;
    cout << "  Everything is logged to: COUT and File." << endl;
    cout << endl;

    test2();

    cout << endl;
    cout << "Test3" << endl;
    cout << "-----" << endl;
    cout << "  Testing Debug system." << endl;
    cout << endl;

    test3();

    return 0;
}

void test1()
{
    // Create stream that would output to COUT and CERR
    //
    MultiStreamBuffer mstream;

    mstream.add(cout.rdbuf());
    mstream.add(cerr.rdbuf());

    Redirector clogRedirector(clog);
    clogRedirector.rdbuf(&mstream);

    clog << "Testing messge that goes into mstream1." << endl;
    clog << "The life is great when everything turns out to be as we expect it to be." << endl;
}

void test2()
{
    // Create stream that outputs information into COUT and FILE
    //
    MultiStreamBuffer mstream;

    std::ofstream fout("test.log");

    mstream.add(cout.rdbuf());
    mstream.add(fout.rdbuf());

    Redirector clogRedirector(clog);
    clogRedirector.rdbuf(&mstream);

    clog << "You mean there are significantly different implementations" << endl;
    clog << "No, I did not mean that" << endl;
    clog << "The iostream library is an object-oriented library." << endl;
}

void test3()
{
    // Create stream that outputs information into COUT and FILE
    //
    Debug debug;
    debug.init("debug.log");

    cout << "Message going into COUT and CLOG." << endl;
    cerr << "CERR and CLOG message." << endl;

    clog << "CLOG only message" << endl;
}
