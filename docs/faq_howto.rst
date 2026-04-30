.. _faq_howto:

FAQ/HOWTO
=========

Retaining coalescence counts
----------------------------

By default, the coalescence counts are temporarily kept in files in the system temp directory, and then removed
prior to ``mrpast process`` completing. If you want to keep these counts for some reason, then set
``MRP_COAL_DIR=<path>`` in your environment when running ``mrpast``.

System temp dir problems
------------------------

If you are having problems with your system temp directory (e.g., on a shared compute cluster), you
can tell ``mrpast`` to use another location by setting ``MRP_TMP_DIR=<path>`` in your environment.