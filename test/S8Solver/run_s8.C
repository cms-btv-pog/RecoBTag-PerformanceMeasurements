#include <memory>

class S8Solver;
std::auto_ptr<S8Solver> s;
// extraOpt = 0 -- NORMAL
// extraOpt = 1 -- mcTrue
// extraOpt = 2 -- ptrel12
// extraOpt = 3 -- ptrel05

// binningChoice = 0 -- {20,30,40,50,60,70,80,90,100,120,140,240}
// binningChoice = 2 --

//        void setbMCscale( const double &value) { _bMCScaleFactor = value; }
//        void setclMCscale( const double &value) { _bMCScaleFactor = value; }
        
void run_s8( const char * data, const char * mc, int extraOpt = 0, double bscale = 1.0, double clscale = 1.0,
                            double pscale = 1.0, double nscale = 1.0, int binningChoice = 1 ) {//, int binCount = 0, Double_t binArray[]) {

    gSystem->Load("libS8Solver.so");

    s.reset(new S8Solver());

    //s->setbMCscale( bscale );
    //s->setclMCscale( clscale );
    //s->setpMCscale( pscale );
    //s->setnMCscale( nscale );


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
    if ( (extraOpt == 1) || (extraOpt == 3) ) {
        s->UseMCTrue(false);
    } else {
        s->UseMCTrue(true);
    }
    if ( extraOpt == 2 ) {
        s->SetPtrelCut( 1.2 );
    }
    if ( extraOpt == 3 ) {
        s->SetPtrelCut( 0.5 );
    }

    switch (binningChoice) {
        case 0: //-- 10GeV 20-100, 20GeV 100-140, overflow @ 140}
            Double_t xbins[12] = { 20,30,40,50,60,70,80,90,100,120,140,240};
            s->SetRebin(12, xbins);
            break;
        case 1:
            Double_t xbins[10] = { 20,30,40,50,60,70,80,100,120,240 };
            s->SetRebin(10, xbins);
            break;
        case 2:
            Double_t xbins[15] = { 20,30,40,50,60,70,80,100,120,160,210,260,320,400,1000 };
            s->SetRebin(14, xbins);
            break;
        default:
            Double_t xbins[12] = { 20,30,40,50,60,70,80,90,100,120,140,240};
            s->SetRebin(12, xbins);
   }
//    s->SetPseudoExperiments( 1000 );
    s->Solve();
    s->PrintData();
    s->Draw();
    s->Save("s8.root");
}
