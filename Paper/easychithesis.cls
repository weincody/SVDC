%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     easychithesis.cls   Version 1.9   Oct  20, 1999
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     Instructions And Hints
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This file provides the LaTeX2e style easychithesis for formating
% graduate thesis at the University of Chicago.  It was written by 
% Bryan Clair and Nathan Dunfield.  For an example its use see the
% file: thesisexample.tex, and the instructions below.  
%
% Comments or problems should be sent to nathand@math.uchicago.edu,
% but I reserve the right to ignore your e-mail, especially after I
% graduate (hopefully June 1999).
%
% Invoking this style:  Start you document with
%
%   \documentclass{easychithesis}
%
% When printing drafts, you may not want your thesis doublespaced.
% If you use instead:
% 
%   \documentclass[singlespace]{easychithesis}
%
% or
%
%   \documentclass[onehalfspace]{easychithesis}
% 
% your thesis will come out single or 3/2 spaced respectively.
% easychithesis also understands all options understood by report.  
% Also the ``double spacing'' provided by this style is not ``true''
% doublespacing as defined by setspace.sty.  Instead, it is the same
% as on the old LaTeX 2.09 thesis style ``chithesis''.  If you want 
% ``true'' double spacing, give the option truedoublespace
%
% Title Page Information Commands
%    Example Usage:
%       \title{Stones And L^2 Invariants}
%       \author{Fred D. Belitnikoff}
%       \department{Mathematics}
%       \division{Physical Sciences}  (Sciences with an s!)
%       \degree{Doctor of Philosophy}
%       \date{July 1998}    (Your grad. date, formatted like the example)
%       \maketitle
%
% \dedication : Use for a dedication, copyright, or epigraph.
%               Produces a page with no number for the text which follows
%               If you want centering, do it yourself with 
%               \begin{center} and \end{center}.
%        Usage:
%           \dedication
%           To blah.
%
% \topmatter : Things like Abstract, Acknowledgements.
%    Usage:
%       \topmatter{Title Of Matter}
%
%       Use like \chapter, once for each page of matter you want.  Don't
%       use numbered sections within these.  Use \section* instead.
%
% The abstract environment has been changed to invoke \topmatter.
%
% \mainmatter : Signals the start of the document's body.
%    Usage:
%       \mainmatter
%
%       Use once before the first numbered \chapter.
%
%  HINTS:
%
%  1.  To get appendices, you don't do anything different from a normal
%     report document.  That means, put the command \appendix before
%     you begin your first appendix, then do each appendix with a
%     \chapter command.  Note that if you have only one appendix, it is
%     customary to leave it unnumbered.  Do this with \chapter*.
%
%  KNOWN PROBLEMS
% 
%  1.  No support for multi-volume theses.
%
%  2.  Doesn't work with older versions of setspace.sty
%
%  3.  Problems with math formulas in chapter headings:
% 
%         a.  Any lowercase letters in the formula are converted to
%         uppercase, e.g. f(x) becomes F(X).   If you really need
%         lowercase math letters in your chapter titles, use the
%         option plainchapterheads (and, if you want, type your
%         chapter titles in ALL CAPS so that the appearance doesn't
%         change).  Note there is no problem 
%         for section or subsection headings in either case.
%
%         b.  Some perfectly reasonable math commands when used in
%         \chapter give the error
%          ``LATEX ERROR: \command  ALLOWED ONLY IN MATH MODE.''
%         The solution to this is to do
%        
%              \newcommand{\mymath}{problem math goes here}
%       
%         and then
%        
%              \chapter{All about \protect\mymath}
%
%         also, the option plainchapterheads will fix this too.
%
%  Improvements over older versions:
%
%  Added in 1.10:  When the table of contents is over 1 page, made so
%  that the page numbers appear on the bottom.
%
%  Added in 1.8:  Increased spacing before Chapter titles.
%
%  Added in 1.7:  Added more space between footnotes.
%
%  Added in 1.6:  Spacing on title page made  made true
%  doublespacing to please thesis office.  This change does not affect
%  the rest of the document.
%
%  Added in 1.4: Name of the bibiograpy is now References, to please the
%  thesis office.  If you want to change it back, put
%
%   \renewcommand{\bibname}{Bibliograpy}   
%
%  in your TeX file prior to \begin{document}


\NeedsTeXFormat{LaTeX2e}[1995/06/01]
\ProvidesClass{easychithesis}[1999/10/20 Chicago Thesis Class For
LaTeX]

\RequirePackage{setspace}[1998/04/8]

\newcommand{\etptsize}{12pt}

\DeclareOption{10pt}{\renewcommand{\etptsize}{10pt}}
\DeclareOption{11pt}{\renewcommand{\etptsize}{11pt}}
\DeclareOption{12pt}{\renewcommand{\etptsize}{12pt}}
\DeclareOption{letterpaper}{\PassOptionsToClass{letterpaper}{report}}

% The following options set the global amount of spacing

\newcommand{\etspacing}{\setstretch{1.37}}
\DeclareOption{truedoublespace}{\renewcommand{\etspacing}{\doublespacing}}
\DeclareOption{doublespace}{\renewcommand{\etspacing}{\setstretch{1.37}}}
\DeclareOption{onehalfspace}{\renewcommand{\etspacing}{\onehalfspacing}}
\DeclareOption{singlespace}{\renewcommand{\etspacing}{\singlespacing}}

% Normally, the chapter heads are all uppercase, but sometimes 
% this cause problems.  The following option fixes changes 
% the chapter heads to small caps instead

\newcommand{\etchapterheadstyle}[1]{\MakeUppercase{#1}}
\DeclareOption{plainchapterheads}{\renewcommand{\etchapterheadstyle}[1]{#1}}

% Pass all other options to report.cls.
\DeclareOption*{%
 \PassOptionsToClass{\CurrentOption}{report}
}

% It would be better to have the following line next, 
% rather than setting the defaults manually, by
% this would break early version of LaTeX2e
%
%\ExecuteOptions{12pt,letterpaper,doublespace}
\ProcessOptions\relax
\LoadClass[\etptsize]{report} 

% Adjust spacing of document

\etspacing

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     Commands Defined For Use By Thesis Author
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% \dedication : Use for a dedication, copyright, or epigraph.
%               Produces a page with no number for the text which follows
%               If you want centering, do it yourself with 
%               \begin{center} and \end{center}.
%
\newcommand{\dedication}{
        \chapter*{}
        \thispagestyle{empty}
}

%
% \topmatter : Things like Abstract, Acknowledgements.
%    Usage:
%       \topmatter{Title Of Matter}
%
%       Use like \chapter, once for each page of matter you want.  Don't
%       use numbered sections within these.  Use \section* instead.
%
\newcommand{\topmatter}[1]{
        \chapter*{#1}
        \addcontentsline{toc}{chapter}{#1}
}

%
% redefinition of abstract environment
% 

\renewenvironment{abstract}{\topmatter{Abstract}}{}

%
% \mainmatter : Signals the start of the document's body.
%    Usage:
%       \mainmatter
%
%       Use once before the first numbered \chapter.
%
\newcommand{\mainmatter}{
        \clearpage
 % In the initial part of the document, all page numbers should be at
 % the bottom.  After \mainmatter  most of them go to the top.
 % Changing where the page numbers go on pages that aren't the first
 % page of a chapter. 
        \pagestyle{myheadings}
        \markright{}
        \pagenumbering{arabic}
}

%
% Title Page Information Commands
%    Example Usage:
%       \title{Stones And L^2 Invariants}
%       \author{Fred D. Belitnikoff}
%       \department{Mathematics}
%       \division{Physical Sciences}  (Sciences with an s!)
%       \date{July 1998}    (Your grad. date, formatted like the example)
%       
%     
\renewcommand{\title}[1]{\gdef\@title{#1}}
\renewcommand{\author}[1]{\gdef\@author{#1}}
\newcommand{\department}[1]{\gdef\@department{#1}}
\newcommand{\division}[1]{\gdef\@division{#1}}
\newcommand{\degree}[1]{\gdef\@degree{#1}}
\renewcommand{\date}[1]{\gdef\@date{#1}}

%
% \maketitle :  Makes the U. of Chicago dissertation title page.
%
\renewcommand{\maketitle}{
   \begin{titlepage}
   \doublespacing
   \begin{center}
   
      THE UNIVERSITY OF CHICAGO
      \vfill
      \MakeUppercase{\@title}\\
      \vfill
      A DISSERTATION SUBMITTED TO\\
      THE FACULTY OF THE DIVISION OF THE \MakeUppercase{\@division}\\
      IN CANDIDACY FOR THE DEGREE OF\\
      \MakeUppercase{\@degree}\\
      \ \\
      DEPARTMENT OF \MakeUppercase{\@department}\\
      \vfill
      BY\\
      \MakeUppercase{\@author}\\
      \vfill
      CHICAGO, ILLINOIS\\
      \MakeUppercase{\@date}
   \end{center}
   \etspacing
   \end{titlepage}
   \pagenumbering{roman}
   \addtocounter{page}{1}
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     Internal Formatting Commands For This Thesis Style
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Redo commands that make chapter & section headers
%
\renewcommand{\@makechapterhead}[1]{{
        \vspace*{0.3in}%
        \normalfont\large\bfseries
        \begin{center}
        \etchapterheadstyle{\@chapapp{ \thechapter}}  %  CHAPTER 1 or APPENDIX A
       
        \etchapterheadstyle{#1}                       %     TITLE OF CHAPTER
        \end{center}
}}
\renewcommand{\@makeschapterhead}[1]{{
    \vspace*{0.3in}%
        \normalfont\large\bfseries
        \begin{center}
        \etchapterheadstyle{#1}                       %     TITLE OF CHAPTER
        \end{center}
}}
\renewcommand\section{\@startsection {section}{1}{\z@}%
                                {-3.5ex \@plus -1ex \@minus -.2ex}%
                                {2.3ex \@plus.2ex}%
                                {\normalfont\large\bfseries\centering}}
\renewcommand\subsection{\@startsection{subsection}{2}{\z@}%
                                {-3.25ex\@plus -1ex \@minus -.2ex}%
                                {1.5ex \@plus .2ex}%
                                {\normalfont\large\itshape\centering}}
\renewcommand\subsubsection{\@startsection{subsubsection}{3}{\z@}%
                                {-3.25ex\@plus -1ex \@minus -.2ex}%
                                {1.5ex \@plus .2ex}%
                                {\normalfont\large}}
\renewcommand\paragraph{\@startsection{paragraph}{4}{\z@}%
                                {3.25ex \@plus1ex \@minus.2ex}%
                                {-1em}%
                                {\normalfont\normalsize\bfseries}}
\renewcommand\subparagraph{\@startsection{subparagraph}{5}{\parindent}%
                                {3.25ex \@plus1ex \@minus .2ex}%
                                {-1em}%
                                {\normalfont\normalsize\itshape}}

% Change style of printing chapters in TOC to match chapter headings.

\renewcommand*\l@chapter[2]{%
  \ifnum \c@tocdepth >\m@ne
    \addpenalty{-\@highpenalty}%
    \vskip 1.0em \@plus\p@
    \setlength\@tempdima{1.5em}%
    \begingroup
      \parindent \z@ \rightskip \@pnumwidth
      \parfillskip -\@pnumwidth
      \leavevmode 
      \advance\leftskip\@tempdima
      \hskip -\leftskip
      \etchapterheadstyle{#1}\nobreak\hfil \nobreak\hb@xt@\@pnumwidth{\hss #2}\par
      \penalty\@highpenalty
    \endgroup
  \fi}

% Change name of table of contents

\renewcommand{\contentsname}{Table of contents}

% Change name of bibliography

\renewcommand{\bibname}{References}

% Change the margins to fit requirements

\setlength{\oddsidemargin}{0.5in}  % really 1.5in
\setlength{\evensidemargin}{0.5in}  % really 1.5in
\setlength{\textwidth}{6in}
\setlength{\topmargin}{0in}   % really 1in 
\setlength{\headheight}{0.20in}
\setlength{\headsep}{0.20in}
\setlength{\textheight}{8.45in}
%\setlength{\footskip}{0.25in}
\setlength{\footnotesep}{0.6cm}



% Make the table of contents, etc.,  single spaced

\let\oldtableofcontents = \tableofcontents
\renewcommand{\tableofcontents}{\begin{singlespace} \oldtableofcontents \end{singlespace}}

\let\oldlistoftables = \listoftables
\renewcommand{\listoftables}{\begin{singlespace} \oldlistoftables \end{singlespace}}

\let\oldlistoffigures = \listoffigures
\renewcommand{\listoffigures}{\begin{singlespace} \oldlistoffigures \end{singlespace}}

%  End of file
