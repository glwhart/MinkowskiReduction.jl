# Dev notes

Artefacts from AI-assisted development sessions on this package,
preserved as teaching examples. They are **not** authoritative
documentation — the current user-facing docs live at
<https://glwhart.github.io/MinkowskiReduction.jl/>.

## What's here

| File | What it is | When it was written |
|---|---|---|
| [`research.md`](research.md) | A deep reading of the algorithm, the 12 Minkowski conditions, the three distinct sources of non-uniqueness in 3D, and the finite-precision pitfalls. Produced by asking the assistant to read the code carefully and also search the web for the Nguyen–Stehlé paper, the Diátaxis framework, and related background. | Session toward v2.0.0 |
| [`plan.md`](plan.md) | Three small, conservative code changes that fell out of the analysis in `research.md` — scale-aware tolerance, a reversed error message, a stale doctest. All three were accepted and landed in commit `d78b511`. | Same session |
| [`../aihelp2-0-1.md`](../aihelp2-0-1.md) | A chronological summary of the prompts and responses across the entire v1.5.0 → v2.0.1 arc. Higher-level than the two files above. | End of the v2.0.1 session |

## How these were produced (for students and collaborators)

Both files were written by Claude Code (Anthropic's Sonnet 4.6 CLI)
during an interactive session, in response to prompts that explicitly
asked for:

1. **Research**: read the actual code plus external sources, produce
   a write-up. This is the kind of task AI is unusually good at because
   it can cross-reference code behaviour with published theory at scale,
   and because the output is a document the human can critique before
   any code is touched.
2. **Concrete change proposals**: "if there are things to change, write
   them in plan.md; otherwise don't". Two-stage output — analysis first,
   proposals second — lets the human sanity-check the analysis before
   deciding whether to act on it.

The value added by careful follow-up prompting (rather than one-shot
generation) was:

- Asking for clarification ("what do you mean by 'load-bearing'?")
  caught phrasing that was clear to the assistant but opaque to a
  reader, and prompted the expanded explanation that ended up in the
  final `research.md`.
- Splitting the work into small review-and-accept cycles ("make these
  three changes", then later "raise the iteration cap", then later
  "restructure the docs") meant every change stood on its own and could
  be reverted without entanglement.
- Insisting on empirical evidence for design decisions (e.g. *"15 was
  motivated, what motivates 50?"*) flushed out analyses that would
  otherwise have stayed implicit, and produced comments in the source
  that document *why* a constant has the value it has.

A rough analogy: the assistant is a tireless junior collaborator who
will write the full first draft of anything you ask for. The value of
the human's time goes entirely into *choosing what to ask for*,
*reviewing what comes back*, and *deciding when to push back*. The two
files in this folder are what the first two of those three do — the
third is what made this package improve rather than drift.

## Not authoritative

If a statement in `research.md` or `plan.md` appears to disagree with
the current source or documentation site, the source/docs win. These
files are snapshots of a moment; the package evolved past them.
