/**
 * LeptonInJetDiscriminator
 * s8 
 *
 * Created by Samvel Khalatian on Nov 17, 2010
 * Copyright 2010, All rights reserved
 */

#ifndef S8_LEPTON_IN_JET_DISCRIMINATOR
#define S8_LEPTON_IN_JET_DISCRIMINATOR

namespace s8
{
    class Lepton;
    class Jet;

    struct BinGroup
    {
        int    bins;
        double min;
        double max;
    };

    // LeptonInJetDiscriminator will return value depending on study
    //
    class LeptonInJetDiscriminator
    {
        public:
            virtual ~LeptonInJetDiscriminator();

            virtual double operator()(const Lepton *, const Jet *) = 0;
            virtual BinGroup bins() const = 0;
    };

    class LeptonInJetPtRelDiscriminator: public LeptonInJetDiscriminator
    {
        public:
            LeptonInJetPtRelDiscriminator();
            virtual ~LeptonInJetPtRelDiscriminator();

            virtual double operator()(const Lepton *, const Jet *);
            virtual BinGroup bins() const;

        private:
            BinGroup _binGroup;
    };

    class LeptonInJetPtDiscriminator: public LeptonInJetDiscriminator
    {
        public:
            LeptonInJetPtDiscriminator();
            virtual ~LeptonInJetPtDiscriminator();

            virtual double operator()(const Lepton *, const Jet *);
            virtual BinGroup bins() const;

        private:
            BinGroup _binGroup;
    };
}

#endif
