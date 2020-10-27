import pyperclip
import re
import sys
from_clipboard = pyperclip.paste()

varX = re.sub(r'1','99',from_clipboard)
varX = re.sub(r'2','1',varX)
varX = re.sub(r'99','2',varX)

vary = re.sub(r'1','99',from_clipboard)
vary = re.sub(r'2','3',vary)
vary = re.sub(r'99','3',vary)

varY = re.sub(r'1','99',from_clipboard)
varY = re.sub(r'2','4',varY)
varY = re.sub(r'99','4',varY)

varz = re.sub(r'1','99',from_clipboard)
varz = re.sub(r'2','5',varz)
varz = re.sub(r'99','5',varz)

varZ = re.sub(r'1','99',from_clipboard)
varZ = re.sub(r'2','6',varZ)
varZ = re.sub(r'99','6',varZ)

varxy = re.sub(r'1','99',from_clipboard)
varxy = re.sub(r'2','7',varxy)
varxy = re.sub(r'99','7',varxy)

varxY = re.sub(r'1','99',from_clipboard)
varxY = re.sub(r'2','8',varxY)
varxY = re.sub(r'99','8',varxY)

varXy = re.sub(r'1','99',from_clipboard)
varXy = re.sub(r'2','9',varXy)
varXy = re.sub(r'99','9',varXy)

varXY = re.sub(r'1','99',from_clipboard)
varXY = re.sub(r'2','10',varXY)
varXY = re.sub(r'99','10',varXY)

varxz = re.sub(r'1','99',from_clipboard)
varxz = re.sub(r'2','11',varxz)
varxz = re.sub(r'99','11',varxz)

varxZ = re.sub(r'1','99',from_clipboard)
varxZ = re.sub(r'2','12',varxZ)
varxZ = re.sub(r'99','12',varxZ)

varXz = re.sub(r'1','99',from_clipboard)
varXz = re.sub(r'2','13',varXz)
varXz = re.sub(r'99','13',varXz)

varXZ = re.sub(r'1','99',from_clipboard)
varXZ = re.sub(r'2','14',varXZ)
varXZ = re.sub(r'99','14',varXZ)

varyz = re.sub(r'1','99',from_clipboard)
varyz = re.sub(r'2','15',varyz)
varyz = re.sub(r'99','15',varyz)

varyZ = re.sub(r'1','99',from_clipboard)
varyZ = re.sub(r'2','16',varZ)
varyZ = re.sub(r'99','16',varZ)

varYz = re.sub(r'1','99',from_clipboard)
varYz = re.sub(r'2','17',varYz)
varYz = re.sub(r'99','17',varYz)

varYZ = re.sub(r'1','99',from_clipboard)
varYZ = re.sub(r'2','18',varYZ)
varYZ = re.sub(r'99','18',varYZ)


lbmBlock = "    // Autogenerating LBM block...\n    " + varX + "\n    " + vary + "\n    " + varY + "\n    " + varz + "\n    " + varZ + "\n    " + varxy + "\n    " + varxY + "\n    " + varXy + "\n    " + varXY + "\n    " + varxz + "\n    " + varxZ + "\n    " + varXz + "\n    " + varXZ + "\n    " + varyz + "\n    " + varyZ + "\n    " + varYz + "\n    " + varYZ + "    //done...\n    "

pyperclip.copy(lbmBlock)
