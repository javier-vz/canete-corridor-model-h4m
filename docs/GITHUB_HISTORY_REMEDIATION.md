# Existing GitHub history: remediation is a separate operation

The currently published `javier-vz/canete-corridor-model` still contains the previous data
organisation, including the old `data/grilla_con_W_hidrica.csv`. Its documentation explicitly
lists the positive archaeological cell flag `presencia_sitio_nuevo` and the pre-revision data
provenance. It also links to historical tags v1.0.0 and v1.0.1.

**A new clean file on `main` will not erase old commit objects, tags, releases, forks,
clone copies or cached views. Do NOT simply commit a new CSV and claim the site locations are
removed.**

Release remediation plan for the repository owner (with project-data authorization):

1. Preserve a secure, access-controlled backup for scientific reproducibility and documentation.
2. Temporarily prevent additional public redistribution while checking all Git refs, tags and
   release assets. Consider setting the repository private if administrator policy permits.
3. In a clean local mirror, identify paths and contents that reveal site locations, including
   site-distance columns, GIS points and derivatives; inspect all refs and tags.
4. Clean historical refs with an appropriate history-rewrite tool (e.g. `git filter-repo`) **only
   after agreeing on exactly which files/content must be removed**. Delete or replace contaminated
   tags/releases/assets and perform an independently checked replacement of the public history.
5. Coordinate GitHub retention/cache and pull-request refs with GitHub Support as needed.
   Force-pushing alone cannot guarantee that previously copied content disappears.
6. Re-clone anonymously and scan all exposed refs and repository-file URLs. Keep a dated audit.
7. If historic exposures cannot be sufficiently controlled, consider a **new clean repository** and
   update the manuscript URL after agreement; this option changes the existing link.

This staging package does not modify GitHub or its refs. Use of the current GitHub URL in the
manuscript is provisional until the public version has been audited.
