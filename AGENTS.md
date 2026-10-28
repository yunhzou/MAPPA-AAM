# Agent instructions

Before starting project work, clone or reuse [PersonalCodingMoat](https://github.com/yunhzou/PersonalCodingMoat.git) and read its [README.md](README.md) and this file in full. Apply these principles alongside the target project's instructions. Reload them for each new session or after context loss.

## 1. Explanation

Keep explanations and documentation short, plain, and understandable to researchers. When communicating with the user, skip niche details unless they matter so the user can focus on important decisions. Say what is needed, then stop. Don’t repeat yourself or keep selling an idea. Explain more when the user asks or when a design decision needs discussion.

## 2. Coding principles

### 2.1 Coding style

Write short, elegant, generalizable code with clear names. Avoid defensive coding and case-specific patches. Fix the underlying logic so the solution applies across cases.

Avoid abbreviations and unnecessary technical terms in names and comments. Follow the explanation rules in section 1: communicate plainly and efficiently for humans. Comments should explain reasoning that the code cannot show.

Keep code easy to follow. Don’t make it clever or compress it just to reduce line count.

### 2.2 Architectural design

Code structure must coevolve with vibe coding. During a rapid proof of concept, it is okay to prototype quickly, even in a single file with hundreds of lines. Once the user confirms that a feature is useful, refactor it into a coherent, modular architecture.

Discuss the high-level design with the developer. Establish the key concepts, module responsibilities, interfaces, and interactions before focusing on niche functionality. Think about the long-term structure and specify which parts of the code belong in each module.

Keep closely related functionality together, minimize dependencies between unrelated components, and separate reusable core logic from experiment-specific scripts. Modularity means meaningful separation of responsibilities, not simply dividing code into smaller files.

Follow the architecture previously approved by the user. Consider how new functionality fits and whether the design needs to change. Avoid competing structures or abstractions without a concrete need.

Revisit the design with the developer when refactoring reveals conflicting responsibilities, functionality does not fit the core modules, or the necessary architectural changes are still unresolved. Resolve these questions before proceeding with the affected refactoring.

You may provide architectural choices, tradeoffs, and recommendations, but the user makes the final decision. Once agreed, record the design in `design.md` and create diagrams showing module responsibilities, hierarchy, interfaces, and interactions.

## 3. Cleaning

Code must never have fallback options. Duplicate functionality is unacceptable. After refactoring, remove the superseded code and options completely from the repository. Leave no obsolete code behind.

When you observe multiple fragmented pieces of code in the repository, treat this as a clear signal that it needs cleaning and refactoring. Follow principle 2.

For example, if several quickly written scripts exist in the repository and the user has confirmed that their functionality is useful, refactor that functionality into reusable core modules instead of leaving it in standalone scripts. Modularity and reusability must be enforced.

## 4. Version control

When the desired outcome is reached and the user confirms that a feature is good, immediately commit it to Git. Clearly document what experiment was performed, what was accepted, and whether refactoring and validation have been completed.

Version design principles as well as code. When the user explicitly specifies how the design should be structured, record those decisions and their reasoning in `design.md` and commit them.

Version control is not limited to preserving code. It must also support effective information exchange by documenting the reasoning and development process. Proper versioning also supports thorough cleaning.

## 5. Working on a cluster

Development often takes place locally and on a cluster. Whenever development spans multiple machines, always use Git to synchronize the work.

There must be a single source of truth for the project, and communication with the user must consistently refer to it. Each machine must have only one copy of the codebase. Enforce this, and use Git to reconcile changes whenever the copies become unsynchronized.
