#include <memory>

class S8Solver;
std::auto_ptr<S8Solver> s;
void run_s8( const char * data, const char * mc, int useMcTrue = 0) {

    gSystem->Load("libS8Solver.so");

    s.reset(new S8Solver());
    s->SetData( data );
    s->SetCorrData( mc );
    s->SetPtrelMaxCut(-1);
    s->SetAlphaConstant(false);
    s->SetBetaConstant(false);
    s->SetGammaConstant(false);
    s->SetDeltaConstant(false);
    s->SetKappabConstant(false);
    s->SetKappaclConstant(false);
    s->setDoAverageSolution(true);
    //s->SetEtaFits();
//    s->setFirstBin(1);
//    s->setLastBin(6);
    //s->SetSolution(0, 1);
    //s->SetSolution(1, 1);
    //s->SetSolution(2, 1);
    //s->SetSolution(3, 1);
    //s->SetSolution(4, 1);
    //s->SetSolution(5, 1);
    //s->SetSolution(6, 1);
    //s->SetSolution(7, 1);
    //s->SetSolution(8, 1);
    //s->SetSolution(9, 1);
    //s->SetSolution(10, 1);
    if (useMcTrue) {
        s->UseMCTrue(false);
    } else {
        s->UseMCTrue(true);
    }
//    s->SetPseudoExperiments( 1000 );
    s->Solve();
    s->PrintData();
    s->Draw();
    s->Save("s8.root");
}
