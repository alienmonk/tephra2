SUBDIRS     = "tephra2"
USER        := $(shell echo ${USER})
REV         := r$(shell svn info | grep "Revision:" | awk '{print $$2}')
USESCRIPT   = "tephra2-${REV}"
ENVIRONDIR  = "/apps/share64/environ"
BIN         = ../bin

all:
	cd tephra2; make

clean:
	cd tephra2; make clean
	rm -f ${USESCRIPT}

install: all
	sed -e "s;@@REV@@;${REV};g" tephra2_use.tmpl > ${USESCRIPT}
	@if [ "${USER}" = "apps" ] ; \
	then \
		if [ -e ${ENVIRONDIR} ] ; \
		then \
			install --mode 0644 -D ${USESCRIPT} ${ENVIRONDIR}; \
		fi \
	fi

distclean: clean
	@if [ "${USER}" = "apps" ] ; \
	then \
		if [ -e "${ENVIRONDIR}/${USESCRIPT}" ] ; \
		then \
			rm -f ${ENVIRONDIR}/${USESCRIPT}; \
		fi \
	fi
	rm -f ${BIN}/tephra2
