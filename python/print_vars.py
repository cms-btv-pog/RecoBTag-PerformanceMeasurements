from variables_cfi import *

variableList.sort(key=lambda x: x.variable)

for variable in variableList:
  name=str(variable.variable)[12:-2]
  description=str(variable.description)[12:-2]
  requires=str(variable.requires)[13:-2]
  if not requires:
    requires='nothing'
  print ('{0: <45}: {1: <80}, requires {2}'.format(name,description,requires))


