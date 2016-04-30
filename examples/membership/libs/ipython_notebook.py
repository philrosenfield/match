"""
Some tools for the notebooks
"""
from IPython.display import display, Markdown
from IPython.nbconvert.filters.markdown import markdown2latex, markdown2html
from IPython.display import DisplayObject
from IPython.html.widgets import FloatProgress
import time as _time
import sys


class Caption(Markdown):
    """ Make a caption to associate with figures """
    def __init__(self, s, center=False, **kwargs):
        Markdown.__init__(self, s, **kwargs)
        self._center = center

    def _repr_html_(self):
        txt = markdown2html(self.data)
        if self._center:
            return '<center>{0}</center>'.format(txt)
        else:
            return '{0}'.format(txt)

    def _repr_latex_(self):
        txt = markdown2latex(self.data)
        if self._center:
            return '\\begin{center}\n' + txt + '\n\\end{center}'
        else:
            return txt

    def display(self):
        display(self)

    def __str__(self):
        return self._repr_latex_()


class Matrix(object):
    """ Make a caption to associate with figures """
    def __init__(self,s, fmt='%0.4g'):
        self.s = s
        self._fmt = fmt

    def _repr_(self):
        text = r"""\begin{bmatrix}"""

        t = []
        for k in self.s:
            t.append( ' & '.join([self._fmt % v for v in k] ) + r'\\' )
        text += ''.join(t)
        text += r"""\end{bmatrix}"""
        return Markdown(text)

    def _repr_latex_(self):
        text = r"""\begin{bmatrix}"""

        t = []
        for k in self.s:
            t.append( ' & '.join([self._fmt % v for v in k] ) + r'\\' )
        text += ''.join(t)
        text += r"""\end{bmatrix}"""
        return text

    def __str__(self):
        return self._repr_latex_()

    def display(self):
        display(self)


def disp_markdown(*args):
    return display(Markdown(*args))


def load_latex_macros():
    return Caption(open('notebook_macros').read(), center=False)


def add_input_toggle():
    from IPython.display import HTML, display

    r = HTML('''
    <script>

    $( document ).ready(function () {
        IPython.CodeCell.options_default['cm_config']['lineWrapping'] = true;
        IPython.notebook.get_selected_cell()

        IPython.toolbar.add_buttons_group([
                {
                    'label'   : 'toggle all input cells',
                    'icon'    : 'fa-eye-slash',
                    'callback': function(){ $('div.input').slideToggle(); }
                }
            ]);
    });
    </script>
    ''')
    display(r)
    return r


def add_citation_button():
    from IPython.display import HTML, display
    r = HTML("""
    <script>

    function insert_citn() {
        // Build paragraphs of cell type and count

        var entry_box = $('<input type="text"/>');
        var body = $('<div><p> Enter the Bibtex reference to insert </p><form>').append(entry_box)
                    .append('</form></div>');

        // Show a modal dialog with the stats
        IPython.dialog.modal({
            notebook: IPython.notebook,
            keyboard_manager: IPython.notebook.keyboard_manager,
            title: "Bibtex reference insertion",
            body: body,
            open: function() {
                // Submit on pressing enter
                var that = $(this);
                that.find('form').submit(function () {
                    that.find('.btn-primary').first().click();
                    return false;
                });
                entry_box.focus();
            },
            buttons : {
                "Cancel" : {},
                "Insert" : {
                    "class" : "btn-primary",
                    "click" : function() {
                        // Retrieve the selected citation, add to metadata,
                        var citation = entry_box.val();
                        // if (!citation) {return;}
                        var citn_html = '<cite data-cite="' + citation + '">' + citation + '</cite>';
                        var cell = IPython.notebook.get_selected_cell();
                        cell.code_mirror.replaceSelection(citn_html);
                    }
                }
            }
        });
    };

    $( document ).ready(function () {

        IPython.toolbar.add_buttons_group([
                {
                    'label'   : 'insert bibtex reference in markdown',
                    'icon'    : 'fa-graduation-cap', // http://fontawesome.io/icons/
                    'callback': insert_citn,
                }
            ]);
    });

    </script>
    <style>
    cite {
        font-style: normal;
        color: #45749e;
    }
    </style>
    """)
    display(r)
    return r


class PDF(object):
    def __init__(self,url):
        self.url = url

    def _repr_html_(self):
        return '<iframe src=%s></iframe>' % self.url

    def _repr_latex_(self):
        return r'\begin{center} \adjustimage{max size={0.9\linewidth}{0.9\paperheight}}{%s}\end{center}' % self.url


class Table(DisplayObject):
    VDOTS = object()

    def __init__(self, data, headings=None, formats=None, caption=None,
                 label=None, position='h', subtables=1):
        """
        A HTML/LaTeX IPython DisplayObject Table

        `data` should be a 2 dimensional array, indexed by row then column,
        with an optional extra row `headings`.

        A 'row' (i.e., an element of `data`) may also be
        :py:const:`Table.VDOTS`, which produces vertical dots in all columns.

        `formats` may be a string, whose format method will be used for every
        cell; a function, called for every cell; or a mixed array of strings
        and functions which is zipped with each row.
        Headings are not formatted.

        `caption` and `label` add the relevant LaTeX markup, and will go in
        the first row of the HTML copy. `label` will have ``tab:`` prepended
        to it.

        If `subtables` is greater than 1, the table will be split into
        `subtables` parts of approximately equal length, and laid out side
        by side.
        """

        if len(data) == 0:
            raise ValueError("data is empty")

        if label is None != caption is None:
            raise ValueError("specify neither or both of label & caption")

        self.columns = len(data[0])
        if self.columns == 0:
            raise ValueError("no columns")

        if headings and len(headings) != self.columns:
            raise ValueError("bad headings length")

        if isinstance(formats, str):
            formats = [formats.format] * self.columns
        elif callable(formats):
            formats = [formats] * self.columns
        elif formats:
            if len(formats) != self.columns:
                raise ValueError("bad formats length")

            def maybe_string_format(f):
                if isinstance(f, str):
                    return f.format
                else:
                    assert callable(f)
                    return f

            formats = list(map(maybe_string_format, formats))
        else:
            formats = [self._default_format] * self.columns

        for i, row in enumerate(data):
            if row is not self.VDOTS and len(row) != self.columns:
                raise ValueError("bad row length", i)

        self.headings = headings
        self.data = data
        self.formats = formats
        self.caption = caption
        self.label = label
        self.position = position
        self.subtables = subtables

    @staticmethod
    def _default_format(what):
        if isinstance(what, float):
            return "{0:.5f}".format(what)
        else:
            return str(what)

    def _format_rows(self):
        for row in self.data:
            if row is self.VDOTS:
                yield self.VDOTS
            else:
                yield (f(x) for f, x in zip(self.formats, row))

    def _subtables_split(self):
        assert self.subtables > 1

        rows = list(self._format_rows())
        nominal_height = len(rows) // self.subtables
        remainder = len(rows) % self.subtables

        heights = [nominal_height] * self.subtables
        for i in range(remainder):
            heights[i] += 1

        slices = []
        acc = 0
        for l in heights:
            slices.append((acc, acc + l))
            acc += l
        assert slices[-1][1] == len(rows)

        subtables = [rows[a:b] for a, b in slices]
        return subtables

    def _repr_latex_(self):
        strings = []

        strings.append(r"""
        \begin{table}[""" + self.position + r"""]
        \centering
        """)

        if self.label:
            strings.append(r"\caption{" + self.caption + "}")
            strings.append(r"\label{tab:" + self.label + "}")

        if self.subtables > 1:
            subtables = self._subtables_split()
            width = "{:.3f}\linewidth".format(0.95 / self.subtables)

            for i, rows in enumerate(subtables):
                strings.append(r"\begin{{subtable}}[t]{{{0}}}%".format(width))
                strings.append(r"""
                \centering
                \vspace{0pt}
                """)
                self._latex_tabular(strings, rows)
                strings.append(r"\end{subtable}%")
                if i != len(subtables) - 1:
                    strings.append("\hfill%")

        else:
            rows = self._format_rows()
            self._latex_tabular(strings, rows)

        strings.append(r"""
        \end{table}
        """)
        return "\n".join(strings)

    def _latex_tabular(self, strings, rows):
        x = "|".join(["c"] * self.columns)
        strings.append(r"\begin{tabular}{|" + x + "|}")
        strings.append(r"\hline")

        if self.headings:
            latex = " & ".join(str(x) for x in self.headings)
            strings.append(latex + r" \\")
            strings.append(r"\hline")

        for row in rows:
            if row is self.VDOTS:
                row = [r"\vdots"] * self.columns
            latex = " & ".join(row)
            strings.append(latex + r" \\")

        strings.append(r"""
        \hline
        \end{tabular}%""")

    def _repr_html_(self):
        strings = []

        strings.append("""
        <style type="text/css">
        .util_Table td { text-align: center; }
        .util_Table tbody tr, .util_Table tbody td {
            border-bottom: 0;
            border-top: 0;
        }
        .util_Table_subtable {
            float: left;
        }
        </style>
        """)

        if self.label:
            c = self.caption
            l = "<code>[{}]</code>".format(self.label)

            strings.append("""
            <h3>{1} {2}</h3>
            """.format(self.columns, c, l))

        if self.subtables > 1:
            subtables = self._subtables_split()
            # width = 0.95 / self.subtables

            strings.append("<div class='clearfix'>")
            for rows in subtables:
                strings.append("<div class='util_Table_subtable'>")
                self._html_table(strings, rows)
                strings.append("</div>")
            strings.append("</div>")

        else:
            rows = self._format_rows()
            self._html_table(strings, rows)

        return "\n".join(strings)

    def _html_table(self, strings, rows):
        strings.append("<table class='util_Table'>")

        if self.headings:
            strings.append("<thead>")
            strings.append("<tr>")
            headings = map("<th>{0}</th>".format, self.headings)
            strings.append("\n".join(headings))
            strings.append("</tr>")
            strings.append("</thead>")

        strings.append("<tbody>")

        for row in rows:
            if row is self.VDOTS:
                row = ["\u22ee"] * self.columns

            strings.append("<tr>")
            row = map("<td>{0}</td>".format, row)
            strings.append("\n".join(row))
            strings.append("</tr>")

        strings.append("</tbody>")
        strings.append("</table>")

    def __repr__(self):
        if self.headings:
            widths = [len(x) for x in self.headings]
            data = [self.headings]
        else:
            widths = None
            data = []

        # don't forget - self._format_rows() is a generator that yields generators
        for row in self._format_rows():
            if row is self.VDOTS:
                continue

            r = list(row)
            w = [len(x) for x in r]

            if widths is None:
                widths = w
            else:
                widths = [max(a, b) for a, b in zip(widths, w)]

            data.append(list(r))

        strings = []
        if self.label:
            c = self.caption.replace("\n", " ")
            strings.append('Table: {0} ({1})'.format(self.label, c))

        for row in data:
            if row is self.VDOTS:
                strings.append('...')
            else:
                r = [x.ljust(b + 4) for x, b in zip(row, widths)]
                strings.append(''.join(r))

        return '\n'.join(strings)

    def __html__(self):
        return self._repr_html_()


class LatexFigure(object):

    extension = 'pdf'

    def __init__(self, label, caption, fig=None, position="", star=False,
                 options='width=\columnwidth'):
        """
        A LaTeX IPython DisplayObject Figure

        `label` is mandatory, since it also sets the filename. It will
        have ``fig:`` preprended to it.

        `fig` is optional - the current figure (via ``gcf``) will be used
        if it is not set.

        `position` is either the float placement specifier or the subfigure
        vertical position.

        If `subfigure` is set to true, a subfigure with width `width` will
        be created.

        The figure is saved (via ``savefig``) as a PDF file in the current
        directory.

        Displaying the object produces LaTeX (only) to embed the figure.
        A little hacky, but since this is meant for use in the notebook
        it is assumed that the figure is going to be displayed automatically
        in HTML independently.
        """
        if fig is None:
            from matplotlib.pyplot import gcf

            fig = gcf()

        self.label = label
        self.caption = caption
        self.fig = fig
        self.position = position
        self.options = options
        self.star = star

        self.filename = "figure_{0:s}.{1:s}".format(label, self.__class__.extension)

        import pylab as plt

        try:
            plt.savefig(self.filename, bbox_inches='tight')
        except:
            plt.savefig(self.filename)

    def _repr_html_(self):
        # Bit crude. Hide ourselves to the notebook viewer, since we'll
        # have been shown already anyway.
        # Nicer solutions are afaict infeasible.
        return markdown2html('> **Figure (<a name="fig:{label:s}">{label:s}</a>)**: {caption:s}'.format(
            label=self.label, caption=self.caption))

    def _repr_latex_(self, subfigure=None):
        if subfigure:
            environment = "subfigure"
            args = "[{position}]{{{width}}}".format(**subfigure)
        else:
            environment = "figure"
            args = "[{0}]".format(self.position)

        if self.star:
            environment += '*'

        return r"""
        \begin{{{env:s}}}{args:s}
            \centering
            \includegraphics[{options:s}]{{{fname:s}}}
            \caption{{{caption:s}}}
            \label{{fig:{label:s}}}
            \end{{{env:s}}}
            """.format(env=environment, args=args, options=self.options,
                       fname=self.filename, caption=self.caption,
                       label=self.label)

        return r"""
        \begin{""" + environment + "}" + args + r"""
            \centering
            \includegraphics[""" + self.options + """]{""" + self.filename + r"""}
            \caption{""" + self.caption + r"""}
            \label{fig:""" + self.label + r"""}
        \end{""" + environment + r"""}
        """

    def __repr__(self):
        c = self.caption.replace("\n", " ")
        return "Figure: {0} ({1})".format(self.label, c)

    def __html__(self):
        return ""


class LatexSubfigures(object):
    def __init__(self, label, caption, figures, position='h',
                 subfigure_position='b'):
        """
        Displays several :cls:`LatexFigures` as sub-figures, two per row.

        `figures` should be an array of :cls:`LatexFigure` objects, not
        :cls:`matplotlib.Figure` objects.
        """

        self.label = label
        self.caption = caption
        self.figures = figures
        self.position = position
        self.subfigure_position = subfigure_position

    def _repr_html_(self):
        # Bit crude. Hide ourselves to the notebook viewer, since we'll
        # have been shown already anyway.
        # Nicer solutions are afaict infeasible.
        return ""

    def _repr_latex_(self):
        strings = []

        strings.append(r"""
        \begin{figure}[""" + self.position + r"""]
            \centering
        """)

        left = True
        first = True
        opts = {"position": self.subfigure_position, "width": ".5\linewidth"}
        for f in self.figures:
            if left and not first:
                strings.append(r"\vspace{1em}")

            # have to be quite careful about whitespace
            latex = f._repr_latex_(subfigure=opts).strip()

            if left:
                latex += '%'
            else:
                latex += r'\newline'

            first = False
            left = not left

            strings.append(latex)

        strings.append(r"""
            \caption{""" + self.caption + r"""}
            \label{fig:""" + self.label + r"""}
        \end{figure}
        """)

        return "\n".join(strings)

    def __repr__(self):
        c = self.caption.replace("\n", " ")
        strings = ["Figure group: {0} ({1})".format(self.label, c)]
        strings += [repr(x) for x in self.figures]
        return "\n".join(strings)

    def __html__(self):
        return ""


class LatexNumberFormatter(object):
    """
    Format floats in exponent notation using latex markup for the exponent

    e.g., ``$-4.234 \\times 10^{-5}$``

    Usage:

    >>> fmtr = LatexNumberFormatter(sf=4)
    >>> fmtr(-4.234e-5)
    "$-4.234 \\\\times 10^{-5}$"
    """

    def __init__(self, sf=10):
        """Create a callable object that formats numbers"""
        self.sf = sf
        self.s_fmt = "{{:.{0}e}}".format(self.sf)

    def __call__(self, n):
        """Format `n`"""
        n = self.s_fmt.format(n)
        n, e, exp = n.partition("e")
        if e == "e":
            exp = int(exp)
            if not n.startswith("-"):
                n = r"\phantom{-}" + n
            return r"${} \times 10^{{{}}}$".format(n, exp)
        else:
            return "${}$".format(n)
