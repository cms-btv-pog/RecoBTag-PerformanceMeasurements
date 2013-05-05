#ifndef TRIGGERPATHS_H 
#define TRIGGERPATHS_H 

#include <map>
#include <string>
#include <boost/assign.hpp>

enum TriggerPaths { 
 HLT_Jet15U,                          
 HLT_Jet30,                           
 HLT_PFJet40,                         
 HLT_Jet30U,                          
 HLT_Jet60,                           
 HLT_Jet50U,                          
 HLT_Jet80,                           
 HLT_PFJet80,                         
 HLT_Jet70U,                          
 HLT_Jet110,                          
 HLT_Jet100U,                         
 HLT_Jet150,                          
 HLT_PFJet140,                        
 HLT_Jet140U,                         
 HLT_Jet190,                          
 HLT_Jet240,                          
 HLT_PFJet260,                        
 HLT_Jet300,                          
 HLT_PFJet320,                        
 HLT_PFJet400,                        
 HLT_DiJetAve15U,                     
 HLT_DiJetAve30,                      
 HLT_DiPFJetAve40,                    
 HLT_DiJetAve30U,                     
 HLT_DiJetAve60,                      
 HLT_DiPFJetAve80,                    
 HLT_DiJetAve50U,                     
 HLT_DiJetAve80,                      
 HLT_DiPFJetAve140,                   
 HLT_BTagMu_Jet10U,                   
 HLT_BTagMu_Jet20U,                   
 HLT_BTagMu_DiJet20U,                 
 HLT_BTagMu_DiJet20U_Mu5,             
 HLT_BTagMu_DiJet20_Mu5,              
 HLT_BTagMu_DiJet30U,                 
 HLT_BTagMu_DiJet30U_Mu5,             
 HLT_BTagMu_DiJet60_Mu7,              
 HLT_BTagMu_DiJet40_Mu5,              
 HLT_BTagMu_DiJet20_L1FastJet_Mu5,    
 HLT_BTagMu_DiJet80_Mu9,              
 HLT_BTagMu_DiJet70_Mu5,              
 HLT_BTagMu_DiJet70_L1FastJet_Mu5,    
 HLT_BTagMu_DiJet100_Mu9,             
 HLT_BTagMu_DiJet110_Mu5,             
 HLT_BTagMu_DiJet110_L1FastJet_Mu5,   
 HLT_BTagMu_Jet300_L1FastJet_Mu5,     
 HLT_BTagMu_Jet300_Mu5 
} ; 
 
std::map<int, std::string> TriggerPathNames = boost::assign::map_list_of 
 (HLT_Jet15U,                          "HLT_Jet15U*")                          
 (HLT_Jet30,                           "HLT_Jet30*")                         
 (HLT_PFJet40,                         "HLT_PFJet40*")                       
 (HLT_Jet30U,                          "HLT_Jet30U*")                          
 (HLT_Jet60,                           "HLT_Jet60*")                         
 (HLT_Jet50U,                          "HLT_Jet50U*")                          
 (HLT_Jet80,                           "HLT_Jet80*")                         
 (HLT_PFJet80,                         "HLT_PFJet80*")                       
 (HLT_Jet70U,                          "HLT_Jet70U*")                          
 (HLT_Jet110,                          "HLT_Jet110*")                        
 (HLT_Jet100U,                         "HLT_Jet100U*")                         
 (HLT_Jet150,                          "HLT_Jet150*")                        
 (HLT_PFJet140,                        "HLT_PFJet140*")                      
 (HLT_Jet140U,                         "HLT_Jet140U*")                         
 (HLT_Jet190,                          "HLT_Jet190*")                        
 (HLT_Jet240,                          "HLT_Jet240*")                        
 (HLT_PFJet260,                        "HLT_PFJet260*")                      
 (HLT_Jet300,                          "HLT_Jet300*")                        
 (HLT_PFJet320,                        "HLT_PFJet320*")                      
 (HLT_PFJet400,                        "HLT_PFJet400*")                      
 (HLT_DiJetAve15U,                     "HLT_DiJetAve15U*")                     
 (HLT_DiJetAve30,                      "HLT_DiJetAve30*")                    
 (HLT_DiPFJetAve40,                    "HLT_DiPFJetAve40*")                  
 (HLT_DiJetAve30U,                     "HLT_DiJetAve30U*")                     
 (HLT_DiJetAve60,                      "HLT_DiJetAve60*")                    
 (HLT_DiPFJetAve80,                    "HLT_DiPFJetAve80*")                  
 (HLT_DiJetAve50U,                     "HLT_DiJetAve50U*")                     
 (HLT_DiJetAve80,                      "HLT_DiJetAve80*")                    
 (HLT_DiPFJetAve140,                   "HLT_DiPFJetAve140*")                 
 (HLT_DiJetAve50U,                     "HLT_DiJetAve50U*")                     
 (HLT_DiJetAve80,                      "HLT_DiJetAve80*")                    
 (HLT_DiPFJetAve140,                   "HLT_DiPFJetAve140*")                 
 (HLT_BTagMu_Jet10U,                   "HLT_BTagMu_Jet10U*")                   
 (HLT_BTagMu_Jet20U,                   "HLT_BTagMu_Jet20U*")                   
 (HLT_BTagMu_DiJet20U,                 "HLT_BTagMu_DiJet20U*")                 
 (HLT_BTagMu_DiJet20U_Mu5,             "HLT_BTagMu_DiJet20U_Mu5*")             
 (HLT_BTagMu_DiJet20_Mu5,              "HLT_BTagMu_DiJet20_Mu5*")              
 (HLT_BTagMu_DiJet30U,                 "HLT_BTagMu_DiJet30U*")               
 (HLT_BTagMu_DiJet30U_Mu5,             "HLT_BTagMu_DiJet30U_Mu5*")             
 (HLT_BTagMu_DiJet60_Mu7,              "HLT_BTagMu_DiJet60_Mu7*")              
 (HLT_BTagMu_DiJet40_Mu5,              "HLT_BTagMu_DiJet40_Mu5*")              
 (HLT_BTagMu_DiJet20_L1FastJet_Mu5,    "HLT_BTagMu_DiJet20_L1FastJet_Mu5*")    
 (HLT_BTagMu_DiJet80_Mu9,              "HLT_BTagMu_DiJet80_Mu9*")              
 (HLT_BTagMu_DiJet70_Mu5,              "HLT_BTagMu_DiJet70_Mu5*")              
 (HLT_BTagMu_DiJet70_L1FastJet_Mu5,    "HLT_BTagMu_DiJet70_L1FastJet_Mu5*")    
 (HLT_BTagMu_DiJet100_Mu9,             "HLT_BTagMu_DiJet100_Mu9*")           
 (HLT_BTagMu_DiJet110_Mu5,             "HLT_BTagMu_DiJet110_Mu5*")             
 (HLT_BTagMu_DiJet110_L1FastJet_Mu5,   "HLT_BTagMu_DiJet110_L1FastJet_Mu5*")   
 (HLT_BTagMu_Jet300_L1FastJet_Mu5,     "HLT_BTagMu_Jet300_L1FastJet_Mu5*")     
 (HLT_BTagMu_Jet300_Mu5,               "HLT_BTagMu_Jet300_Mu5*") ;                
 
enum TriggerPaths_TTBar {
HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL, 
HLT_Mu17_Mu8,                                                                         
HLT_Mu17_TkMu8,                                      
HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL,  
HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL   
} ; 

std::map<int, std::string> TriggerPathNames_TTBar = boost::assign::map_list_of 
  (HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL, "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL*") 
  (HLT_Mu17_Mu8,                                                                         "HLT_Mu17_Mu8*") 
  (HLT_Mu17_TkMu8,                                                                       "HLT_Mu17_TkMu8*")                                  
  (HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL,                                     "HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL*")
  (HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL,                                     "HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL*") ; 

#endif 
