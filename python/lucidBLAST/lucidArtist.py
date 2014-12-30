__author__ = 'Hunter Cameron'

# TODO Add a method to print results instead of displaying.
# TODO Add the option to plot more than one subject at once. Perhaps using subplots 2x2?
        # TODO The subplots should have the ability to be switched out. Perhaps a subplot list in one pannel?
# TODO The alignment that is selected should be a different color than the other alignments


# TODO Make plotting window appear maximized from the start.

import matplotlib.pyplot as plt
import lucidSubject
import lucidAlignment




class MultiPlotScaffold(object):
    """
    Class to handle the plotting of lucidBLAST scaffolds.

    Call the plot_scaffolds method to plot
    """

    def __init__(self, scaffolds, height=1000, width=1000, num_panels=3):
        self.scaffolds = scaffolds
        self.height = height
        self.width = width
        self.plt_lists = {}
        self.dragged_artist = None  # initialize variable

        # if the user requested more panels than subjects, trim the num panels to the num subjects
        if len(scaffolds) < num_panels:
            self.num_panels = len(scaffolds)
        elif num_panels > 5:
            self.num_panels = 5
        else:
            self.num_panels = num_panels

        # set up figure (entire plotting window)
        self.fig = plt.figure("lucidBLAST Results", tight_layout=True)
        # add mouse listeners
        self.fig.canvas.mpl_connect('pick_event', self.on_pick)
        self.fig.canvas.mpl_connect('button_release_event', self.on_release)

        # make navigation pannel
        nav = plt.subplot2grid((self.num_panels, 5), (0, 0), rowspan=self.num_panels)
        # remove axes
        nav.axes.get_xaxis().set_visible(False)
        nav.axes.get_yaxis().set_visible(False)
        nav.set_title("Subject List")
        nav.label_dict = {}
        self.nav = nav

        # make plotting panels
        self.axes = []
        for num in range(self.num_panels):

            ax = plt.subplot2grid((self.num_panels, 5), (num, 1), colspan=4)
            ax.axes.get_yaxis().set_visible(False)
            ax.tick_params(axis='x', which='both', bottom='off', top='off')     # remove axis ticks for the x axis
            self.axes.append(ax)

    def plot_scaffolds(self):
        """
        Driver method to plot the instance scaffolds.
        :return:
        """
        self.make_plotting_lists()

        # populate navigation pane and plotting panes
        text_spacer = .05       # TODO place the text boxes by some sort of absolute scale, maybe make y axis 1:1000?  <- Won't work but something like that. Some large number to keep it consistent at different window sizes.
        for index, name in enumerate(sorted(self.plt_lists)):
            textbox = self.nav.text(.05, .99 - (text_spacer * index), name, label=self.scaffolds[name][0], picker=5, transform=self.nav.transAxes, verticalalignment="top")
            textbox.set_clip_on(True)
            self.nav.label_dict[str(self.scaffolds[name][0])] = self.scaffolds[name][0]

            if index < self.num_panels:     # plot if there is an avail pane
                self.plot_list(self.axes[index], name)

        plt.show()

    def plot_list(self, ax, name):
        """
        Plots the scaffold contained in self.scaffolds[name] to the Axes object ax
        :param ax: Pyplot Axes object
        :param name: name string
        :return:
        """

        # set the axis to render the plot in the middle
        length = self.scaffolds[name][0].length
        xmin = int(length * -.1)
        xmax = int(length + (length * .1))
        ymin = -10
        ymax = 10
        ax.axis([xmin, xmax, ymin, ymax])

        # reset the annotation
        ax.annotation = ax.text(0, -1, '', verticalalignment="top", family="monospace")
        ax.annotation.set_clip_on(True)


        ax.label_dict = {}
        for element in self.plt_lists[name]:
            label, x, y, col = element
            name_x, name_y, obj = label
            # Add title
            if isinstance(obj, lucidSubject.Subject):
                #ax.text(0, ymax - 1, "Subject: {}".format(obj.name), verticalalignment="top")
                ax.set_title(obj.name)
            ax.label_dict[str(obj)] = obj     # labels are converted to text. Needs to way to get back to object.

            ax.plot(x, y, col, linewidth=5, label=obj, picker=5)

    def on_pick(self, event):
        """
        Simple method to store the artist that was clicked for use in the on_release method.
        :param event:
        :return:
        """
        artist = event.artist
        self.dragged_artist = artist

    def on_release(self, event):
        """
        Fired when the mouse button is released. Processes what should be done as a result of the click.
        :param event:
        :return:
        """

        if self.dragged_artist is None:     # if click began in the void...
            pass
        elif event.inaxes is None:          # if the click was released in the void...
            pass
        elif event.inaxes == self.nav:      # no changing the navigation pane!
            pass

        elif event.inaxes == self.dragged_artist.get_axes():
            self.annotation_handler()
        else:
            self.switch_handler(event.inaxes)

        # do I need to reset this? Yes, because mouse releases will be fired even when pick events aren't first.
        self.dragged_artist = None

    def switch_handler(self, new_axes):
        """
        Plots a new scaffold in the Axes new_axes. Scaffold that will be plotted is stored in self.dragged_artist.
        :param new_axis: A matplotlib pyplot Axes object.
        :return:
        """

        # currently only want to switch if the label was a subject.
        label = self.dragged_artist.get_label()
        obj_ref = self.dragged_artist.get_axes().label_dict[label]

        if isinstance(obj_ref, lucidSubject.Subject):
            # make a new subplot
            new_axes.cla()
            self.plot_list(new_axes, obj_ref.name)

            self.fig.canvas.draw()

    def annotation_handler(self):
        """
        Gets and outputs the annotation for a clicked object.
        :return:
        """

        artist = self.dragged_artist
        ax = artist.get_axes()          # This claims err in editor but it's fine.
        label = artist.get_label()      # Same

        obj_ref = ax.label_dict[label]      # TODO some sort of handler for key not in dict error? Technically don't need becuase key should always be there, right?
        string = ''
        if isinstance(obj_ref, lucidSubject.Subject):
            string = "\n".join(["Subject Sequence",
                                "Name:   {}".format(obj_ref.name),
                                "Length: {}".format(obj_ref.length), ])

        elif isinstance(obj_ref, lucidAlignment.Alignment):
            string = "\n".join(["Alignment to Subject: {}".format(obj_ref.subject.name),
                                "Align Length:     {}".format(obj_ref.length),
                                "Query Name:       {}".format(obj_ref.query.name),
                                "Query Length:     {}".format(obj_ref.query.length),
                                "Percent Identity: {}".format(obj_ref.perc_identity),
                                "Evalue:           {}".format(obj_ref.evalue),
                                "Query Start:      {}".format(obj_ref.qry_start),
                                "Query End:        {}".format(obj_ref.qry_end),
                                "Subj Start:       {}".format(obj_ref.subj_start),
                                "Subj End:         {}".format(obj_ref.subj_end), ])

        else:
            string = 'Label not found. Likely an internal error with naming conventions'

        ax.annotation.set_text(string)
        self.fig.canvas.draw()
        print(string)

    def make_plotting_lists(self):
        """
        Converts a scaffold into a list of line segments to plot.
        :return:
        """
        for name, scaffold in self.scaffolds.items():

            plt_list = []
            for tier in range(len(scaffold)):
                if tier == 0:       # plot subject line
                    ln_tuple = ((0, 0, scaffold[0]), (0, scaffold[0].length), (tier, tier), "black")
                    plt_list.append(ln_tuple)
                else:               # plot alignments
                    for alignment in scaffold[tier]:
                        if alignment.reverse:
                            color = "red"
                        else:
                            color = "blue"
                        #color = "red" if alignment.reverse == True else color = "blue"     # conditional doesn't work. why?
                        ln_tuple = ((alignment.subj_start, tier, alignment), (alignment.subj_start, alignment.subj_end), (tier, tier), color)
                        plt_list.append(ln_tuple)

            self.plt_lists[name] = plt_list
