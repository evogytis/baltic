# Contributing to baltic

baltic is an open source/science tool and very open to contributions that fix and/or expand its utility. While phylogenetc data manipulation is one aspect of baltic it is understood that visualising such data will be the main application and as such the design philosophy of baltic follows closely that of matplotlib - functionality in baltic itself should be basic and generally applicable with more involved visualisations combining elements thereof. The most sought-after contributions are therefore bug fixes, code to manipulate trees and the occasional novel visualisation technique. 

Note that baltic follows and enforces the [Contributor Covenant Code of Conduct](https://www.contributor-covenant.org/).

## Contributions

If you'd like to contribute a bug fix feel free to start a pull request outlining the problem and its solution.
For adding code for a new tree manipulation function or for a novel visualisation please open an issue first and describe the goal and suggested approach first.
Even if suggested manipulations/visualisations don't make it into baltic's code base (_e.g._ if the application isn't broad) there will be room to share it under examples.

#### Release process  
Notes for internal use:  
1. Make a new github tag following [semantic versioning](https://medium.com/the-non-traditional-developer/semantic-versioning-for-dummies-45c7fe04a1f8)  
2. Edit the [version number](https://github.com/evogytis/baltic/blob/master/setup.py#L11) and [download url](https://github.com/evogytis/baltic/blob/master/setup.py#L14) in `setup.py`  
3. Push to pypi

```
pip install twine  
python3 setup.py sdist 
twine upload dist/*
```  
