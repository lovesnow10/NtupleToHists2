#include "TTreeFormulaContainer.hpp"
#include <TIterator.h>

using namespace std;

TTreeFormulaContainer::TTreeFormulaContainer(bool setOwner) {
  tlist = new TList();
  tlist->SetOwner(setOwner);
}

bool TTreeFormulaContainer::Notify(void) {

  bool success = true;

  TIter ti(tlist);
  TObject *temp;

  while ((temp = ti.Next())) {
    ((TTreeFormula *)temp)->UpdateFormulaLeaves();
    success = temp->Notify() && success;
  }

  return success;
}

bool TTreeFormulaContainer::AddFormula(TTreeFormula *ttf) {

  int n = tlist->GetEntries();
  tlist->Add(ttf);

  return (tlist->GetEntries() == n + 1);
}

bool TTreeFormulaContainer::UnsetNotify(TTreeFormula *ttf) {

  Int_t n = tlist->GetEntries();
  while (tlist->Remove(ttf))
    ; // Removes all copies of ttf

  return (tlist->GetEntries() < n);
}

TTreeFormulaContainer::~TTreeFormulaContainer() { delete tlist; }
