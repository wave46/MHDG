# Build identity embedded in each solver executable.
# Regression builds provide these values explicitly. Direct make builds derive
# stable values from the repository and remain buildable outside a Git checkout.
MHDG_REPOSITORY_ROOT ?= $(abspath $(CURDIR)/..)
MHDG_GIT_COMMIT ?= $(shell \
	cd "$(MHDG_REPOSITORY_ROOT)" 2>/dev/null \
	&& git rev-parse HEAD 2>/dev/null \
	|| printf unknown)
MHDG_GIT_DIRTY ?= $(shell \
	if cd "$(MHDG_REPOSITORY_ROOT)" 2>/dev/null \
		&& git rev-parse \
		--is-inside-work-tree >/dev/null 2>&1; then \
	  if test -n "$$(git status \
			--porcelain --untracked-files=all)"; then \
	    printf true; \
	  else \
	    printf false; \
	  fi; \
	else \
	  printf false; \
	fi)
MHDG_BUILD_ID ?= manual-$(MHDG_GIT_COMMIT)$(if $(filter true,$(MHDG_GIT_DIRTY)),-dirty)

BUILD_PROVENANCE_SOURCE := build_provenance.F
BUILD_PROVENANCE_OBJECT := build_provenance.o
BUILD_PROVENANCE_DIRTY := $(if $(filter true,$(MHDG_GIT_DIRTY)),.TRUE.,.FALSE.)

.PHONY: force-build-provenance
force-build-provenance:

$(BUILD_PROVENANCE_SOURCE): force-build-provenance
	@{ \
	  printf '%s\n' 'MODULE build_provenance'; \
	  printf '%s\n' '  IMPLICIT NONE'; \
	  printf '%s\n' '  PRIVATE'; \
	  printf '%s\n' '  PUBLIC :: solver_git_commit, solver_git_dirty, solver_build_id'; \
	  printf '%s\n' "  CHARACTER(LEN=*), PARAMETER :: solver_git_commit = '$(MHDG_GIT_COMMIT)'"; \
	  printf '%s\n' "  LOGICAL, PARAMETER :: solver_git_dirty = $(BUILD_PROVENANCE_DIRTY)"; \
	  printf '%s\n' "  CHARACTER(LEN=*), PARAMETER :: solver_build_id = '$(MHDG_BUILD_ID)'"; \
	  printf '%s\n' 'END MODULE build_provenance'; \
	} > $@.tmp
	@if $(CMP) -s $@.tmp $@; then \
	  $(RM) $@.tmp; \
	else \
	  $(CP) $@.tmp $@; \
	  $(RM) $@.tmp; \
	fi

$(BUILD_PROVENANCE_OBJECT): $(BUILD_PROVENANCE_SOURCE)
	$(FC) -c $(FCFLAGS) $(INC) $< -o $@
