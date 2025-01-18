<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
<xsl:output method="text"/>
 <xsl:strip-space elements="*"/>

<xsl:template match="Test/Results/Measurement/Value/text()">
  <xsl:value-of select="concat(.,'&#xA;')"/>
</xsl:template>

<xsl:template match="text()"/>

</xsl:stylesheet>