/*
  A simple class to store TTreeFormula
  Modified based on https://root.cern.ch/phpBB3/viewtopic.php?t=15549
*/
#ifndef TTREEFORMULAGCONTAINER_HPP_
#define TTREEFORMULAGCONTAINER_HPP_

#include <TObject.h>
#include <TTreeFormula.h>
#include <TList.h>

using namespace std;

class TTreeFormulaContainer : public TObject {

public:
  // Constructor: If setOwner==kTRUE, deletes member TTreeFormulas on delete.
  TTreeFormulaContainer(bool setOwner = true);

  // Destructor: Self-explanatory.
  virtual ~TTreeFormulaContainer();

  // Notify: Calls UpdateFormulaLeaves on all member TTreeFormulas
  // Also calls Notify for each member TTreeFormula
  // Returns true if every member Notify returns true
  bool Notify(void);

  // SetNotify: Adds passed TTreeFormula to group
  // Returns true if TTreeFormula successfully added.
  bool AddFormula(TTreeFormula *);

  // UnsetNotify: Removes passed TTreeFormula from group.
  // Returns true if TTreeFormula successfully removed.
  // Returns false if it was not a member.
  bool UnsetNotify(TTreeFormula *);

  long GetEntries() { return tlist->GetEntries(); };
  TTreeFormula *GetFormula(int ientry) {
    TTreeFormula *f = (TTreeFormula *)tlist->At(ientry);
    return f;
  }

private:
  // Pay no attention to the TList behind the curtain.
  TList *tlist;
};

#endif /* TTREEFORMULAGROUP_H_ */
