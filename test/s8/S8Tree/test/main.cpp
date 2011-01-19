#include <iostream>
#include <string>
#include <vector>
#include <utility>

#include <boost/functional/hash.hpp>

#include <TRandom3.h>
#include <TClass.h>

#include "interface/S8Event.h"
#include "interface/S8Trigger.h"
#include "interface/S8TriggerCenter.h"

using std::cout;
using std::endl;
using std::make_pair;
using std::string;
using std::vector;

using s8::TriggerCenter;

typedef vector<string> TriggerNames;
typedef boost::hash<string> MyHash;

std::ostream &operator<<(std::ostream &out, const TriggerCenter &triggerCenter);
void merge(TObject *, const TObject *);

void fillTriggers(TriggerCenter &triggerCenter, TRandom &);

int main(int argc, char **argv)
{
    TRandom3 randomizer;

    TriggerCenter triggerCenter;
    fillTriggers(triggerCenter, randomizer);

    cout << triggerCenter << endl;

    TriggerCenter triggerCenter2;
    fillTriggers(triggerCenter2, randomizer);

    cout << triggerCenter2 << endl;

    TriggerCenter triggerCenterMerge;
    cout << "Merge" << endl;
    cout << triggerCenterMerge << endl;

    merge(&triggerCenterMerge, &triggerCenter);
    cout << "Merge" << endl;
    cout << triggerCenterMerge << endl;

    merge(&triggerCenterMerge, &triggerCenter2);
    cout << "Merge" << endl;
    cout << triggerCenterMerge << endl;

    return 1;
}

std::ostream &operator<<(std::ostream &out, const TriggerCenter &triggerCenter)
{
    const TriggerCenter::TriggerMap &triggers = triggerCenter.triggers();

    for(TriggerCenter::TriggerMap::const_iterator trigger = triggers.begin();
        triggers.end() != trigger;
        ++trigger)
    {
        out << trigger->first << ": " << trigger->second << endl;
    }

    return out;
}

void merge(TObject *to, const TObject *from)
{
    if (from->IsA()->InheritsFrom(TriggerCenter::Class()) &&
        to->IsA()->InheritsFrom(TriggerCenter::Class()))
    {
        cout << "Safe to merge" << endl;

        dynamic_cast<TriggerCenter *>(to)->merge(*dynamic_cast<const TriggerCenter *>(from));
    }
    else
        cout << " NOT safe to merge" << endl;
}

void fillTriggers(TriggerCenter &triggerCenter, TRandom &random)
{
    TriggerNames triggerNames;
    triggerNames.push_back("HLT1");
    triggerNames.push_back("HLT2");
    triggerNames.push_back("HLT3");
    triggerNames.push_back("Adium");
    triggerNames.push_back("Test");
    triggerNames.push_back("AlfredApp");
    triggerNames.push_back("AMC");
    triggerNames.push_back("Peter");
    triggerNames.push_back("Jackson");
    triggerNames.push_back("Battery");
    triggerNames.push_back("Crazy");
    triggerNames.push_back("Another one");
    triggerNames.push_back("BOOST ROCKS");
    triggerNames.push_back("Peace of Cake");
    triggerNames.push_back("What else");
    triggerNames.push_back("Honestly");
    triggerNames.push_back("I don't know what to add");
    triggerNames.push_back("Different Space");
    triggerNames.push_back("Omni Trio");
    triggerNames.push_back("Safari");
    triggerNames.push_back("Apple");
    triggerNames.push_back("Steve");
    triggerNames.push_back("Jobs");
    triggerNames.push_back("Chrome");
    triggerNames.push_back("Music");
    triggerNames.push_back("LTJ Buckem");
    triggerNames.push_back("Paper");
    triggerNames.push_back("Piper");
    triggerNames.push_back("Window");
    triggerNames.push_back("Mouse");
    triggerNames.push_back("Phone");
    triggerNames.push_back("Scotch");
    triggerNames.push_back("Drum");
    triggerNames.push_back("Growl");

    TriggerCenter::TriggerMap &triggers = triggerCenter.triggers();

    const int max = random.Uniform(1, 10);
    cout << "Fill " << max << " entries" << endl;
    for(int i = 0; max > i; ++i)
    {
        const int id = random.Uniform(0, triggerNames.size());

        MyHash hash;
        std::size_t h = hash(triggerNames[id]);
        triggers[h] = triggerNames[id];
    }
}
