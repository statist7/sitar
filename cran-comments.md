## R CMD check results
There were no ERRORs or WARNINGs.

There was a NOTE of a missing URL (which I didn't see):

   Found the following (possibly) invalid URLs:
     URL: ...
       From: man/deren.Rd
             man/ob_convertr.Rd
             man/who0607.Rd
       Status: 404
       Message: Not Found

I have now updated the URL in the three files and there are no NOTEs.

## Downstream dependencies
I have run revdepcheck for downstream dependencies.
MonoInc the one relevant package passed.
