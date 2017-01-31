<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
<xsl:output method="xml" indent="yes"/>

<xsl:template match="/Catch">
  <testsuites>
      <xsl:attribute name="name"><xsl:value-of select="@name"/></xsl:attribute>
      <xsl:attribute name="tests"><xsl:value-of select="OverallResults/@successes"/></xsl:attribute>
      <xsl:attribute name="failures"><xsl:value-of select="OverallResults/@failures"/></xsl:attribute>
      <xsl:apply-templates/>
  </testsuites>
</xsl:template>

<xsl:template match="TestCase">
  <testsuite>
    <xsl:attribute name="name"><xsl:value-of select="@name"/></xsl:attribute>
    <xsl:apply-templates/>
  </testsuite>
</xsl:template>

<xsl:template match="Section">
  <testcase>
      <xsl:attribute name="name"><xsl:value-of select="@name"/></xsl:attribute>
      <xsl:attribute name="assertions"><xsl:value-of select="OverallResults/@successes"/></xsl:attribute>
      <xsl:apply-templates/>
  </testcase>
</xsl:template>

<xsl:template match="Info">
  <system-out><xsl:value-of select="text()"/></system-out>
</xsl:template>

<xsl:template match="Expression">
   <error>
   <xsl:attribute name="type"><xsl:value-of select="@type"/></xsl:attribute>
   <xsl:attribute name="message">
       <xsl:value-of select="@filename"/>:<xsl:value-of select="@line"/>: 
       <xsl:value-of select="Expanded/text()"/>
   </xsl:attribute>
   </error>
</xsl:template>

</xsl:stylesheet>
