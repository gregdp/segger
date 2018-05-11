from  chimera.extension import EMO, manager

# -----------------------------------------------------------------------------
#
class Segment_Map_EMO ( EMO ):

  def name(self):
    return 'Segment Map'
  def description(self):
    return self.categoryDescriptions()['Volume Data']
  def categories(self):
    return self.categoryDescriptions().keys()
  def categoryDescriptions(self):
    # since we want to use specialized descriptions for certain categories...
    return {
      'Volume Data': 'Segment volume data to identify components',
    }
  def icon(self):
    return self.path('volseg.png')
  def activate(self):
    # self.module('volumedialog').show_volume_dialog()
    d = self.module('segment_dialog').show_volume_segmentation_dialog()
    return None

# -----------------------------------------------------------------------------
#
class Fit_Segments_EMO ( EMO ):

  def name(self):
    return 'Fit to Segments'
  def description(self):
    return self.categoryDescriptions()['Volume Data']
  def categories(self):
    return self.categoryDescriptions().keys()
  def categoryDescriptions(self):
    # since we want to use specialized descriptions for certain categories...
    return {
      'Volume Data': 'Fit structures into segmented regions',
    }
  def icon(self):
    return None
  def activate(self):
    # self.module('volumedialog').show_volume_dialog()
    d = self.module('fit_dialog').show_fit_segments_dialog()
    return None

# -----------------------------------------------------------------------------
# Register dialogs and menu entry.
#
manager.registerExtension ( Fit_Segments_EMO ( __file__ ) )
manager.registerExtension ( Segment_Map_EMO ( __file__ ) )

# -----------------------------------------------------------------------------
# Register segmentation file reader.
#
def open_seg(path):
  from Segger import segment_dialog as sd
  d = sd.show_volume_segmentation_dialog()
  d.OpenSegFiles([(path, 'Segmentation')])
  return []

import chimera
fi = chimera.fileInfo
fi.register('Segger segmentation', open_seg, ['.seg'], ['segger'],
            category = fi.GENERIC3D)

# -----------------------------------------------------------------------------
# Register segment command.
#
def segment(cmdname, args):
  from Segger import segcmd
  segcmd.segment_command(cmdname, args)

import Midas.midas_text
Midas.midas_text.addCommand('segment', segment, None, help = True)
