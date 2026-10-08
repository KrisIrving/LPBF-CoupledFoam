# Project continuity

Before continuing development, read `docs/project-status.zh-CN.md`,
`docs/research-plan.zh-CN.md`, and the latest entries in
`docs/test-register.zh-CN.md` and `docs/development-log.md`.

The user runs compilation and CFD tests. Do not invoke WSL or run simulations
on the development machine without new user authorization. Reading supplied
feedback archives is allowed. Develop, inspect statically, document and push
within the user's authorized project scope.

Use milestone development packages and integrated test matrices. Do not
continually add isolated uniform-case micro-tests. M1-A is closed as a limited
audit, not universal physical validation. Next work is M1-B as described in the
research plan. Q01 remains unresolved; more outer correctors did not fix it.

For each feedback archive, record its tested commit, archive SHA-256, test
scope, outcome, numerical evidence, limitations and changed hypotheses in
`docs/test-results.json` and `docs/test-register.zh-CN.md`. Update project
status and development log; update the plan when scope or decisions change.
Do not overwrite historical failures with later pass results. Missing or
unevaluated checks are not passes. Do not claim unexecuted runtime validation.

Keep raw reports, machine details, credentials and private manuscripts out of
the public repository. Curated result values and project-authored documentation
may be committed. Do not alter numerical acceptance thresholds after observing
results merely to make a test pass.
