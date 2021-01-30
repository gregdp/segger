from  chimera.extension import EMO, manager

# -----------------------------------------------------------------------------
#
class Segger_EMO ( EMO ):

  def name(self):
    return 'Segger'
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
class SegFit_EMO ( EMO ):

  def name(self):
    return 'SegFit'
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
    return self.path('fitseg.png')
  def activate(self):
    # self.module('volumedialog').show_volume_dialog()
    d = self.module('fit_dialog').show_fit_segments_dialog()
    return None

# -----------------------------------------------------------------------------
#
class MapQ_EMO ( EMO ):

  def name(self):
    return 'MapQ from Segger'
  def description(self):
    return self.categoryDescriptions()['Volume Data']
  def categories(self):
    return self.categoryDescriptions().keys()
  def categoryDescriptions(self):
    # since we want to use specialized descriptions for certain categories...
    return {
      'Volume Data': 'Evaluate map & model',
    }
  def icon(self):
    return self.path('mapq.png')
  def activate(self):
    # self.module('volumedialog').show_volume_dialog()
    d = self.module('mapq').show_dialog()
    return None


# -----------------------------------------------------------------------------
# Register dialogs and menu entry.
#
manager.registerExtension ( SegFit_EMO ( __file__ ) )
manager.registerExtension ( Segger_EMO ( __file__ ) )
#manager.registerExtension ( MapQ_EMO ( __file__ ) )

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
