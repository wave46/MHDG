"""Extract the production stabilization and prescribed-trace assembly for a focused check."""
from pathlib import Path
import sys

jacobian = Path(sys.argv[1]).read_text()
boundary = Path(sys.argv[2]).read_text()
face = jacobian.split('    SUBROUTINE assemblyExtFacesContribution(', 1)[1]
stabilization = face.split('      ! Assembly stabilization terms', 1)[1].split(
    '      !************* End stabilization terms', 1)[0]
Path('dirichlet_stabilization.F').write_text(stabilization)
start = boundary.index('  SUBROUTINE assembly_dirichletstr_bc(')
end = boundary.index('    ! Compute Q^T^', start)
Path('dirichlet_prescribed.F').write_text(
    boundary[start:end] + '  END SUBROUTINE assembly_dirichletstr_bc\n')
