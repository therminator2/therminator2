#ifndef HGEANT_EVENT_SAVER
    #define HGEANT_EVENT_SAVER

    #include "AbstractEventSaver.h"

    class HGeantEventSaver : public AbstractEventSaver 
    {
        public:
            HGeantEventSaver();
            virtual void Save(Event *, Model *, int);
        private:
            int GetGeantID(int therminatorPID);
            std::map<int,int> PIDmap;
    };

#endif
