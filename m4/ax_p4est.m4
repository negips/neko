AC_DEFUN([AX_P4EST],[
    # p4est paths
    AC_ARG_WITH([p4est],
    AS_HELP_STRING([--with-p4est=DIR],
    [Compile with support for p4est]),
    [
    if test -d "$withval"; then
       ac_p4est_path="$withval";
       ac_p4est_libdir="$ac_p4est_path/lib"
       ac_p4est_incdir="$ac_p4est_path/include"
    fi		
    ],[with_p4est=no])

    # libsc path if different from p4est
    AC_ARG_WITH([libsc],
    AS_HELP_STRING([--with-libsc=DIR],
    [Path to libsc if different from p4est path]),
    [
    if test -d "$withval"; then
       ac_libsc_path="$withval";
       ac_libsc_libdir="$ac_libsc_path/lib"
       ac_libsc_incdir="$ac_libsc_path/include"
    fi
    ],[with_libsc=no])

    if test -d "$ac_p4est_libdir"; then
       P4EST_LDFLAGS="-L$ac_p4est_libdir"
    else
	   with_p4est=no
    fi

    if test -d "$ac_p4est_incdir"; then
       P4EST_CPPFLAGS="-I$ac_p4est_incdir"
    else
       with_p4est=no
    fi

    if test "x${with_libsc}" != xno; then
       if test -d "$ac_libsc_libdir"; then
          LIBSC_LDFLAGS="-L$ac_libsc_libdir"
       fi
       if test -d "$ac_libsc_incdir"; then
          LIBSC_CPPFLAGS="-I$ac_libsc_incdir"
       fi
    fi

    if test "x${with_p4est}" != xno; then
       CPPFLAGS_SAVED="$CPPFLAGS"
       LDFLAGS_SAVED="$LDFLAGS"
       CPPFLAGS="$P4EST_CPPFLAGS $CPPFLAGS"
       LDFLAGS="$P4EST_LDFLAGS $LDFLAGS"
       if test "x${with_libsc}" != xno; then
          CPPFLAGS="$LIBSC_CPPFLAGS $CPPFLAGS"
          LDFLAGS="$LIBSC_LDFLAGS $LDFLAGS"
       fi
       export CPPFLAGS
       export LDFLAGS
       
       have_p4est=yes

       AC_LANG(C)
       # test p4est
       AC_CHECK_HEADER([p4est.h],[have_p4est_h=yes],[have_p4est_h=no])
       if test x"${have_p4est_h}" = xno; then
          CPPFLAGS="$CPPFLAGS_SAVED"
          have_p4est=no
          AC_MSG_WARN([p4est.h not found])
       fi

       AC_CHECK_LIB(p4est, p4est_init,
                    [have_p4est_l=yes;P4EST_LIBS="-lp4est"],
                    [have_p4est_l=no],[-lp4est])
       if test x"${have_p4est_l}" = xno; then
          LDFLAGS="$LDFLAGS_SAVED"
          P4EST_LIBS=""
          have_p4est=no
          AC_MSG_WARN([p4est_init in libp4est not found])
       fi

       # test libsc
       AC_CHECK_HEADER([sc.h],[have_libsc_h=yes],[have_libsc_h=no])
       if test x"${have_libsc_h}" = xno; then
          CPPFLAGS="$CPPFLAGS_SAVED"
          have_p4est=no
          AC_MSG_WARN([sc.h not found])
       fi

       AC_CHECK_LIB(sc, sc_init,
                    [have_libsc=yes;P4EST_LIBS="$P4EST_LIBS -lsc"],
                    [have_libsc=no],[-lsc])
       if test x"${have_libsc}" = xno; then
          LDFLAGS="$LDFLAGS_SAVED"
          P4EST_LIBS=""
          have_p4est=no
          AC_MSG_WARN([sc_init in libsc not found])
       fi
       AC_LANG(Fortran)
       AC_SUBST(P4EST_LIBS)

       # pre-processing flag
       if test x"${have_p4est}" = xyes; then
          AC_DEFINE(HAVE_P4EST,1,[Define if you have the p4est and sc libraries.])
       fi
    fi
])

