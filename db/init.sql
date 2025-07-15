-- phpMyAdmin SQL Dump
-- version 5.2.1
-- https://www.phpmyadmin.net/
--
-- Host: 127.0.0.1:28901
-- Generation Time: Jul 14, 2025 at 07:40 AM
-- Server version: 10.11.5-MariaDB-log
-- PHP Version: 8.2.4

SET SQL_MODE = "NO_AUTO_VALUE_ON_ZERO";
START TRANSACTION;
SET time_zone = "+00:00";


/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!40101 SET NAMES utf8mb4 */;

--
-- Database: `dtuser_DASH`
--
DROP DATABASE IF EXISTS `dtuser_DASH`;
CREATE DATABASE IF NOT EXISTS `dtuser_DASH` DEFAULT CHARACTER SET utf8mb3 COLLATE utf8mb3_unicode_ci;
USE `dtuser_DASH`;

DELIMITER $$
--
-- Procedures
--
DROP PROCEDURE IF EXISTS `Del_Basecall`$$
CREATE DEFINER=`dtuser`@`172.%` PROCEDURE `Del_Basecall` (IN `lids` TEXT)   BEGIN
	DELETE FROM Basecall
	WHERE LID IN (lids);
END$$

DROP PROCEDURE IF EXISTS `Del_BasecallPeaks`$$
CREATE DEFINER=`dtuser`@`172.%` PROCEDURE `Del_BasecallPeaks` (IN `lids` TEXT)   BEGIN
	DELETE FROM Basecall_Peaks
	WHERE LID IN (lids);
END$$

DROP PROCEDURE IF EXISTS `Del_CentroidDistance`$$
CREATE DEFINER=`dtuser`@`172.%` PROCEDURE `Del_CentroidDistance` (IN `lids` TEXT, IN `clusters` TEXT, IN `method` VARCHAR(100))   BEGIN
	DELETE FROM Centroid_Distance
	WHERE FIND_IN_SET(Centroid_Distance.LID, lids) > 0 AND FIND_IN_SET(Centroid_Distance.cluster, clusters) > 0 AND Centroid_Distance.`method`=`method`;
END$$

DROP PROCEDURE IF EXISTS `Del_CentroidRegions`$$
CREATE DEFINER=`dtuser`@`172.%` PROCEDURE `Del_CentroidRegions` (IN `lids` TEXT, IN `method` VARCHAR(100))   BEGIN
	DELETE FROM Centroid_Regions
	WHERE LID IN (lids) AND Centroid_Regions.method=`method`;
END$$

DROP PROCEDURE IF EXISTS `Del_Centroids`$$
CREATE DEFINER=`dtuser`@`172.%` PROCEDURE `Del_Centroids` (IN `seq_ids` TEXT, IN `method` VARCHAR(100))   BEGIN
	DELETE FROM Centroids
	WHERE FIND_IN_SET(Centroids.LID, seq_ids) > 0 AND Centroids.`method`=method;
END$$

DROP PROCEDURE IF EXISTS `Del_CentroidSecondary`$$
CREATE DEFINER=`dtuser`@`172.%` PROCEDURE `Del_CentroidSecondary` (IN `lids` TEXT, IN `clusters` TEXT, IN `method` VARCHAR(100))   BEGIN
	DELETE FROM Centroid_Secondary
	WHERE LID IN (lids) AND cluster IN (clusters) AND Centroid_Secondary.method=`method`;
END$$

DROP PROCEDURE IF EXISTS `Del_CentroidSecondaryIntrx`$$
CREATE DEFINER=`dtuser`@`172.%` PROCEDURE `Del_CentroidSecondaryIntrx` (IN `lid1` INT, IN `lid2` INT, IN `cluster` INT, IN `cluster2` INT, IN `method` VARCHAR(100))   BEGIN
	SELECT @ssid:=SSID FROM Centroid_Secondary_Interactions ssi
	WHERE ssi.LID1=lid1 AND ssi.LID2=lid2 AND ssi.cluster1=cluster AND ssi.cluster2=cluster2 AND ssi.method=method;

	CASE WHEN @ssid IS NOT NULL THEN
		BEGIN
			DELETE FROM Centroid_Secondary_Interactions
			WHERE SSID IN (@ssid);
		
			DELETE FROM Centroid_Probabilities
			WHERE SSID IN (@ssid);
		
			DELETE FROM Centroid_BPP
			WHERE SSID IN (@ssid);
		END;
	ELSE BEGIN END;
	END CASE;
END$$

DROP PROCEDURE IF EXISTS `Del_Clusters`$$
CREATE DEFINER=`dtuser`@`172.%` PROCEDURE `Del_Clusters` (IN `seq_ids` TEXT, IN `method` VARCHAR(100))   BEGIN
	DELETE FROM Clusters
	WHERE FIND_IN_SET(Clusters.LID, seq_ids) > 0 AND Clusters.`method`=method;
END$$

DROP PROCEDURE IF EXISTS `Del_Gmm`$$
CREATE DEFINER=`dtuser`@`172.%` PROCEDURE `Del_Gmm` (IN `lids` TEXT)   BEGIN
	DELETE FROM Gmm
	WHERE LID IN (lids);
END$$

DROP PROCEDURE IF EXISTS `Del_Lof`$$
CREATE DEFINER=`dtuser`@`172.%` PROCEDURE `Del_Lof` (IN `lids` TEXT)   BEGIN
	DELETE FROM Lof
	WHERE LID IN (lids);
END$$

DROP PROCEDURE IF EXISTS `Del_Peaks`$$
CREATE DEFINER=`dtuser`@`172.%` PROCEDURE `Del_Peaks` (IN `seq_ids` TEXT)   BEGIN
	DELETE FROM Peaks
	WHERE FIND_IN_SET(Peaks.LID, seq_ids) > 0;
END$$

DROP PROCEDURE IF EXISTS `Del_ReactivityFull`$$
CREATE DEFINER=`dtuser`@`172.%` PROCEDURE `Del_ReactivityFull` (IN `lids` TEXT)   BEGIN
	DELETE FROM Reactivity_full WHERE LID IN (lids);
END$$

DROP PROCEDURE IF EXISTS `Del_ReadDepth`$$
CREATE DEFINER=`dtuser`@`172.%` PROCEDURE `Del_ReadDepth` (IN `lids` TEXT)   BEGIN
	DELETE FROM Read_depth WHERE LID IN (lids);
END$$

DROP PROCEDURE IF EXISTS `Del_ReadDepthFull`$$
CREATE DEFINER=`dtuser`@`172.%` PROCEDURE `Del_ReadDepthFull` (IN `lids` TEXT)   BEGIN
	DELETE FROM Read_depth_full WHERE LID IN (lids);
END$$

DROP PROCEDURE IF EXISTS `Del_ReadDepthFullRead`$$
CREATE DEFINER=`dtuser`@`172.%` PROCEDURE `Del_ReadDepthFullRead` (IN `lid` TEXT, IN `read` TEXT)   BEGIN
	DELETE FROM Read_depth_full
	WHERE LID IN (`lid`) AND read_index IN (`read`);
END$$

DROP PROCEDURE IF EXISTS `Del_Signal`$$
CREATE DEFINER=`dtuser`@`172.%` PROCEDURE `Del_Signal` (IN `lids` TEXT)   BEGIN
	DELETE FROM Modified
	WHERE LID  IN (lids);

	DELETE FROM Unmodified
	WHERE LID  IN (lids);
END$$

DROP PROCEDURE IF EXISTS `Del_StructureBpp`$$
CREATE DEFINER=`dtuser`@`172.%` PROCEDURE `Del_StructureBpp` (IN `ssid` UUID)   BEGIN
	DELETE FROM Structure_Probabilities
	WHERE SSID=ssid;

	DELETE FROM Structure_BPP
	WHERE SSID=ssid;
END$$

DROP PROCEDURE IF EXISTS `Del_StructureSecondaryIntrx`$$
CREATE DEFINER=`dtuser`@`172.%` PROCEDURE `Del_StructureSecondaryIntrx` (IN `lid1` INT, IN `lid2` INT, IN `read_index` INT)   BEGIN
	SELECT @ssid:=SSID FROM Structure_Secondary_Interactions ssi
	WHERE ssi.LID=lid1 AND ssi.LID2=lid2 AND ssi.read_index=read_index;

	CASE WHEN @ssid IS NOT NULL THEN
		BEGIN
			DELETE FROM Structure_Secondary_Interactions
			WHERE SSID=@ssid;
		
			DELETE FROM Structure_Probabilities
			WHERE SSID=@ssid;
		
			DELETE FROM Structure_BPP
			WHERE SSID=@ssid;
		END;
	ELSE BEGIN END;
	END CASE;
END$$

DROP PROCEDURE IF EXISTS `Gmm_Control`$$
CREATE DEFINER=`dtuser`@`172.%` PROCEDURE `Gmm_Control` (IN `lids` TEXT, IN `types` VARCHAR(100))   BEGIN
	SELECT g.LID, g.`position`, g.Predict as 'gmm', sc.predict as 'control'
	FROM Gmm g
	INNER JOIN Structure_Library sl ON sl.LID=g.LID 
	INNER JOIN Structure_Control sc ON sc.SID=sl.SID AND sc.`position`=g.`position`
	WHERE g.LID IN (lids) AND FIND_IN_SET(sc.type,`types`);
END$$

DROP PROCEDURE IF EXISTS `Gmm_Structure`$$
CREATE DEFINER=`dtuser`@`172.%` PROCEDURE `Gmm_Structure` (IN `lids` TEXT)   BEGIN
	SELECT g.LID, g.`position`, g.Predict as 'gmm', s.base_type, s.metric
	FROM Gmm g
	INNER JOIN Structure_Library sl ON sl.LID=g.LID 
	INNER JOIN Structure s ON s.SID=sl.SID AND s.`position`=g.`position`
	WHERE g.LID IN (lids);
END$$

DROP PROCEDURE IF EXISTS `Ins_Library`$$
CREATE DEFINER=`dtuser`@`172.%` PROCEDURE `Ins_Library` (IN `vcontig` VARCHAR(100), IN `vsequence` TEXT, IN `vsequence_name` VARCHAR(250), IN `vsequence_len` INT, IN `vtemp` FLOAT, IN `vtype1` VARCHAR(50), IN `vtype2` VARCHAR(50), IN `vcomplex` INT, IN `vis_modified` INT, IN `vis_putative` BINARY, IN `vdatetime` DATETIME, IN `vrun` INT, OUT `lid` TEXT)   BEGIN
	INSERT INTO Library (contig, sequence, sequence_name, sequence_len, temp, type1, type2, complex, is_modified, is_putative, timestamp, run)
	VALUES (vcontig, vsequence, vsequence_name, vsequence_len, vtemp, vtype1, vtype2, vcomplex, vis_modified, vis_putative,vdatetime, vrun);
	SELECT @lid := LAST_INSERT_ID();
END$$

DROP PROCEDURE IF EXISTS `Ins_Mfe`$$
CREATE DEFINER=`dtuser`@`172.%` PROCEDURE `Ins_Mfe` (IN `id` INT, IN `lid` INT, IN `mfe` FLOAT, IN `timestamp` DATETIME)   BEGIN
	UPDATE Structure_Secondary SET Structure_Secondary.mfe= mfe, Structure_Secondary.timestamp =`timestamp`
	WHERE Structure_Secondary.ID=id AND Structure_Secondary.LID=lid;
END$$

DROP PROCEDURE IF EXISTS `Ins_StructureLibrary`$$
CREATE DEFINER=`dtuser`@`172.%` PROCEDURE `Ins_StructureLibrary` (IN `lid` INT, IN `contigs` VARCHAR(50), OUT `msid` INT)   BEGIN
	-- insert structure id, same for mod and unmodifed, by get max current sid
	SELECT MAX(SID) INTO @msid FROM Structure_Library WHERE contig=contigs; 
	IF @msid IS NULL 
	THEN
		SELECT MAX(SID)+1 INTO @msid FROM Structure_Library;
	END IF;
	-- insert into structure library joiner table
	INSERT INTO Structure_Library (LID, SID, contig)
	VALUES (lid, @msid,contigs);	
	SELECT @msid;
END$$

DROP PROCEDURE IF EXISTS `SelectBC`$$
CREATE DEFINER=`dtuser`@`172.%` PROCEDURE `SelectBC` (IN `lids` TEXT)   BEGIN
	SELECT d.LID, l.contig, `position`, insertion, mismatch, deletion, quality, basecall_reactivity, aligned_reads, d.`sequence`, l.run, l.temp, l.complex, l.type1, l.type2
	FROM Basecall d
	INNER JOIN Library l ON l.ID=d.LID
	WHERE l.ID IN (lids)
	ORDER BY contig,`position`;
END$$

DROP PROCEDURE IF EXISTS `SelectBCUnmod`$$
CREATE DEFINER=`dtuser`@`172.%` PROCEDURE `SelectBCUnmod` (IN `lids` TEXT)   BEGIN
	SELECT d.LID, d.contig, `position`, AVG(insertion) as insertion, AVG(mismatch) as mismatch, AVG(deletion) as deletion, AVG(quality) as quality, AVG(basecall_reactivity) as basecall_reactivity, AVG(aligned_reads) as aligned_reads, d.sequence, l.temp, l.complex, l.type1, l.type2
	FROM Basecall d
	INNER JOIN Library l ON l.ID=d.LID
	WHERE l.ID IN (lids)
	GROUP BY d.LID, d.contig, `position`, d.`sequence`
	ORDER BY d.LID, d.contig, `position`;
END$$

DROP PROCEDURE IF EXISTS `SelectCentroids`$$
CREATE DEFINER=`dtuser`@`172.%` PROCEDURE `SelectCentroids` (IN `lids` TEXT)   BEGIN
	SELECT c.LID, c.contig, position, cluster, centroid as 'reactivity', method, l.sequence
	FROM Centroids c 
	INNER JOIN Library l ON l.ID=c.LID
	WHERE c.LID IN (lids)
	ORDER BY c.LID, c.contig, position, cluster, method;
END$$

DROP PROCEDURE IF EXISTS `SelectClusters`$$
CREATE DEFINER=`dtuser`@`172.%` PROCEDURE `SelectClusters` (IN `contig` VARCHAR(255), IN `temp` INT, IN `type1` VARCHAR(255), IN `type2` VARCHAR(255), IN `complex` INT, IN `method` VARCHAR(100))   BEGIN
	SELECT contig, position, cluster, centroid, method
	FROM Centroids c 
	INNER JOIN Library l ON l.LID=c.LID
	WHERE l.contig=contig AND l.temp=temp AND l.type1=type1 AND l.type2=type2 AND l.complex=complex AND c.method=method
	ORDER BY contig, position, cluster;
END$$

DROP PROCEDURE IF EXISTS `SelectLibrary`$$
CREATE DEFINER=`dtuser`@`172.%` PROCEDURE `SelectLibrary` (IN `contig` VARCHAR(255), IN `temp` INT, IN `type1` VARCHAR(255), IN `type2` VARCHAR(255), IN `complex` INT)   BEGIN
	SELECT *
	FROM Library l
	WHERE l.contig=contig AND l.temp=temp AND l.type1=type1 AND l.type2=type2 AND l.complex=complex;
END$$

DROP PROCEDURE IF EXISTS `SelectLibraryById`$$
CREATE DEFINER=`dtuser`@`172.%` PROCEDURE `SelectLibraryById` (IN `lids` TEXT)   BEGIN
	SELECT l.*, ss.secondary, ss.MFE
	FROM Library l
	LEFT JOIN Structure_Secondary ss ON ss.LID=l.ID
	WHERE l.ID IN (lids);
END$$

DROP PROCEDURE IF EXISTS `SelectLibraryLids`$$
CREATE DEFINER=`dtuser`@`172.%` PROCEDURE `SelectLibraryLids` (IN `contig` VARCHAR(255), IN `temp` INT, IN `type1` VARCHAR(255), IN `type2` VARCHAR(255), IN `complex` INT)   BEGIN
	SELECT DISTINCT ID
	FROM Library l
	WHERE l.contig=contig AND l.temp=temp AND l.type1=type1 AND l.type2=type2 AND l.complex=complex;
END$$

DROP PROCEDURE IF EXISTS `SelectMaxStructureProbabilities`$$
CREATE DEFINER=`dtuser`@`172.%` PROCEDURE `SelectMaxStructureProbabilities` (IN `ssid` BIGINT)   BEGIN
	SELECT LID, base as 'position', read_index, MAX(max) as 'base_pair_prob' 
	FROM (SELECT LID, base1 as 'base', read_index, MAX(probability) as 'max' 
	FROM  Structure_Secondary_Interactions ssi
	INNER JOIN Structure_BPP sb ON sb.SSID=ssi.SSID
	INNER JOIN Structure_Probabilities sp  ON sp.SSID=ssi.SSID
	WHERE ssi.SSID=ssid AND ssi.type='MFE' AND type_interaction IN ('A', 'AA')
	GROUP BY base1, type_interaction
	UNION 
	SELECT LID, base2 as 'base', read_index, MAX(probability)
	FROM  Structure_Secondary_Interactions ssi
	INNER JOIN Structure_BPP sb ON sb.SSID=ssi.SSID
	INNER JOIN Structure_Probabilities sp  ON sp.SSID=ssi.SSID
	WHERE ssi.SSID=ssid AND ssi.type='MFE' AND type_interaction IN ('A', 'AA')
	GROUP BY base1, type_interaction) AS bpp
	GROUP BY bpp.base;
END$$

DROP PROCEDURE IF EXISTS `SelectMod`$$
CREATE DEFINER=`dtuser`@`172.%` PROCEDURE `SelectMod` (IN `contig` VARCHAR(255), IN `temp` INT, IN `type1` VARCHAR(255), IN `type2` VARCHAR(255), IN `complex` INT)   BEGIN
	SELECT d.LID, d.contig, read_index, `position`, event_level_mean, event_length, l.run, l.temp, l.complex, l.type1, l.type2
	FROM Library l
	INNER JOIN Modified d ON d.LID=l.ID
	WHERE l.contig=contig AND l.temp=temp AND l.type1=type1 AND l.type2=type2 AND l.complex=complex
	ORDER BY contig, read_index, `position`, run, temp;
END$$

DROP PROCEDURE IF EXISTS `SelectPeaks`$$
CREATE DEFINER=`dtuser`@`172.%` PROCEDURE `SelectPeaks` (IN `mod_lids` TEXT, IN `unmod_lids` TEXT)   BEGIN
	DROP TEMPORARY TABLE IF EXISTS t1;

	CREATE TEMPORARY TABLE t1 AS
		SELECT contig, `position`, AVG(event_length) AS 'event_length_dmso', AVG(event_level_mean) AS 'event_level_mean_dmso'
		FROM Unmodified d
		WHERE d.LID IN (unmod_lids)
		GROUP BY contig, `position`;
	
	SELECT a.LID, a.contig, a.position, a.read_index, (ABS(event_length) - ABS(event_length_dmso)) as 'delta_dwell', (ABS(event_level_mean) - ABS(event_level_mean_dmso)) as 'delta_signal'   
	FROM Modified a
	INNER JOIN t1 ON t1.contig=a.contig AND t1.position=a.position
	WHERE a.LID IN (mod_lids)
	ORDER BY a.LID, a.contig, a.read_index, a.position;
END$$

DROP PROCEDURE IF EXISTS `SelectPredict`$$
CREATE DEFINER=`dtuser`@`172.%` PROCEDURE `SelectPredict` (IN `lids` TEXT)   BEGIN
	DROP TEMPORARY TABLE IF EXISTS t1;
	DROP TEMPORARY TABLE IF EXISTS s1;
	DROP TEMPORARY TABLE IF EXISTS s2;
	
	CREATE TEMPORARY TABLE t1 AS
		SELECT DISTINCT read_index FROM Peaks s
		WHERE s.LID IN (lids);
	
	CREATE TEMPORARY TABLE s1 AS
		SELECT DISTINCT LID, s.contig, position, sequence FROM `Structure` s 
		LEFT JOIN Structure_Library sl ON sl.SID=s.SID
		WHERE LID IN (lids);
	
	CREATE TEMPORARY TABLE s2
		(UNIQUE struct (LID, contig, read_index, position)) AS
		SELECT LID, contig, read_index, position, sequence FROM t1
		CROSS JOIN s1;

	SELECT * FROM (SELECT St.LID, St.contig, St.read_index, St.position, St.sequence as 'Sequence', b.predict as 'Predict_BC',  p.predict_signal as 'Predict_Signal', p.predict_dwell as 'Predict_Dwell', ld.predict_dwell as 'Predict_Lofd', ld.predict_signal as 'Predict_Lofs', g.predict as 'Predict_Gmm'
	FROM s2 as St
	LEFT JOIN Basecall_Peaks b ON b.LID=St.LID AND b.`position`=St.`position`
	LEFT JOIN Gmm g ON g.LID=St.LID AND g.`position`=St.`position`
	LEFT JOIN Peaks p ON p.LID=St.LID AND p.read_index=St.read_index AND p.`position`=St.`position`
	LEFT JOIN Lof ld ON ld.LID=St.LID  AND ld.read_index=St.read_index  AND ld.`position`=St.`position` 
	WHERE St.LID IN (lids)) f
	ORDER BY LID, read_index, `position`;
END$$

DROP PROCEDURE IF EXISTS `SelectPredictControl`$$
CREATE DEFINER=`dtuser`@`172.%` PROCEDURE `SelectPredictControl` (IN `lids` TEXT)   BEGIN
	DROP TEMPORARY TABLE IF EXISTS t1;
	DROP TEMPORARY TABLE IF EXISTS s1;
	DROP TEMPORARY TABLE IF EXISTS s2;
	
	CREATE TEMPORARY TABLE t1 AS
		SELECT DISTINCT read_index FROM Peaks s
		WHERE s.LID IN (lids);
	
	CREATE TEMPORARY TABLE s1 AS
		SELECT DISTINCT LID, contig, position, sequence FROM `Structure` s 
		WHERE s.LID IN (lids);
	
	CREATE TEMPORARY TABLE s2
		(UNIQUE struct (LID, contig, read_index, position)) AS
		SELECT LID, contig, read_index, position, sequence FROM t1
		CROSS JOIN s1;

	SELECT * FROM (SELECT St.LID, St.contig, St.read_index, St.position, St.sequence as 'Sequence', b.predict as 'Predict_BC',  p.predict_signal as 'Predict_Signal', p.predict_dwell as 'Predict_Dwell', ld.predict_dwell as 'Predict_Lofd', ld.predict_signal as 'Predict_Lofs', g.predict as 'Predict_Gmm'
	FROM s2 as St
	LEFT JOIN Basecall_Peaks b ON b.LID=St.LID AND b.`position`=St.`position`
	LEFT JOIN Gmm g ON g.LID=St.LID AND g.`position`=St.`position`
	LEFT JOIN Peaks p ON p.LID=St.LID AND p.read_index=St.read_index AND p.`position`=St.`position`
	LEFT JOIN Lof ld ON ld.LID=St.LID  AND ld.read_index=St.read_index  AND ld.`position`=St.`position` 
	WHERE St.LID IN (lids)) f
	ORDER BY contig, read_index, `position`;
END$$

DROP PROCEDURE IF EXISTS `SelectReadDepthFull`$$
CREATE DEFINER=`dtuser`@`172.%` PROCEDURE `SelectReadDepthFull` (IN `lids` TEXT)   BEGIN
	DROP TEMPORARY TABLE IF EXISTS T2;
	CREATE TEMPORARY TABLE IF NOT EXISTS T2 AS 
	(SELECT read_index, COUNT(position) as 'cp'
	FROM Modified m 
	WHERE LID IN (lids)
	GROUP BY read_index
	UNION 
	SELECT read_index, COUNT(position) as 'cp'
	FROM Unmodified u 
	WHERE LID IN (lids)
	GROUP BY read_index);

	SELECT AVG(cp) INTO @avgcp FROM T2;
	
	DROP TEMPORARY TABLE IF EXISTS LR;
	CREATE TEMPORARY TABLE IF NOT EXISTS LR AS 
	(SELECT read_index FROM T2
	WHERE T2.cp >= @avgcp);
	
	SELECT St.contig, CONCAT(St.read_index,St.LID) as read_id, St.read_index, St.position, St.sequence as 'Sequence', St.Reactivity_Score as 'Reactivity', St.Predict, St.LID
	FROM Read_depth_full as St
	INNER JOIN Library l ON l.ID=St.LID
	INNER JOIN LR ON LR.read_index=St.read_index
	WHERE St.LID IN (lids)
	ORDER BY St.contig, read_index, `position`, l.run;
END$$

DROP PROCEDURE IF EXISTS `SelectStructure`$$
CREATE DEFINER=`dtuser`@`172.%` PROCEDURE `SelectStructure` (IN `lids` TEXT)   BEGIN
	SELECT s.ID, sl.LID, s.position, s.contig, s.sequence, s.base_type
	FROM Structure s
	INNER JOIN Structure_Library sl ON sl.SID=s.SID
	WHERE sl.LID IN (lids);
END$$

DROP PROCEDURE IF EXISTS `SelectUnmod`$$
CREATE DEFINER=`dtuser`@`172.%` PROCEDURE `SelectUnmod` (IN `contig` VARCHAR(255), IN `temp` INT, IN `type1` VARCHAR(255), IN `type2` VARCHAR(255), IN `complex` INT)   BEGIN
	SELECT d.LID, d.contig, read_index, `position`, event_level_mean, event_length, l.run, l.temp, l.complex, l.type1, l.type2
	FROM Library l
	INNER JOIN Unmodified d ON d.LID=l.ID
	WHERE l.contig=contig AND l.temp=temp AND l.type1=type1 AND l.type2=type2 AND l.complex=complex
	ORDER BY contig, read_index, `position`, run, temp;
END$$

DROP PROCEDURE IF EXISTS `Select_BcControl`$$
CREATE DEFINER=`dtuser`@`172.%` PROCEDURE `Select_BcControl` (IN `lids` TEXT, IN `types` VARCHAR(100))   BEGIN
	SELECT g.LID, g.`position`, g.predict, sc.predict as 'control'
	FROM Basecall_Peaks g
	INNER JOIN Structure_Library sl ON sl.LID=g.LID 
	INNER JOIN Structure_Control sc ON sc.SID=sl.SID AND sc.`position`=g.`position`
	WHERE g.LID IN (lids) AND FIND_IN_SET(sc.type,`types`);
END$$

DROP PROCEDURE IF EXISTS `Select_BcStructure`$$
CREATE DEFINER=`dtuser`@`172.%` PROCEDURE `Select_BcStructure` (IN `lids` TEXT)   BEGIN
	SELECT g.LID, g.`position`, g.predict, s.base_type, s.metric
	FROM Basecall_Peaks g
	INNER JOIN Structure_Library sl ON sl.LID=g.LID 
	INNER JOIN Structure s ON s.SID=sl.SID AND s.`position`=g.`position`
	WHERE g.LID IN (lids);
END$$

DROP PROCEDURE IF EXISTS `Select_Centroids`$$
CREATE DEFINER=`dtuser`@`172.%` PROCEDURE `Select_Centroids` (IN `lids` TEXT)   BEGIN
	SELECT c.LID, c.contig, c.position, c.cluster, c.centroid, l.sequence,c.`method`, cs.secondary, ss.secondary as 'control_secondary'
	FROM Centroids c
	LEFT JOIN Centroid_Secondary cs ON cs.LID=c.LID AND cs.cluster=c.cluster AND cs.method=c.method
	LEFT JOIN Structure_Secondary ss ON ss.LID=c.LID
	INNER JOIN Library l ON l.ID=c.LID
	WHERE c.LID IN (lids) AND cs.type="MFE" OR cs.type IS NULL
	ORDER BY LID, c.contig, position, cluster, c.method;
END$$

DROP PROCEDURE IF EXISTS `Select_DwellControl`$$
CREATE DEFINER=`dtuser`@`172.%` PROCEDURE `Select_DwellControl` (IN `lids` TEXT, IN `types` VARCHAR(100))   BEGIN
	SELECT g.LID, g.read_index, g.`position`, g.predict_dwell as 'predict', sc.predict as 'control'
	FROM Peaks g
	INNER JOIN Structure_Library sl ON sl.LID=g.LID 
	INNER JOIN Structure_Control sc ON sc.SID=sl.SID AND sc.`position`=g.`position`
	WHERE g.LID IN (lids) AND FIND_IN_SET(sc.type,`types`);
END$$

DROP PROCEDURE IF EXISTS `Select_DwellStructure`$$
CREATE DEFINER=`dtuser`@`172.%` PROCEDURE `Select_DwellStructure` (IN `lids` TEXT)   BEGIN
	SELECT g.LID, g.read_index, g.`position`, g.predict_dwell as 'predict', s.base_type, s.metric
	FROM Peaks g
	INNER JOIN Structure_Library sl ON sl.LID=g.LID 
	INNER JOIN Structure s ON s.SID=sl.SID AND s.`position`=g.`position`
	WHERE g.LID IN (lids);
END$$

DROP PROCEDURE IF EXISTS `Select_LibraryFull`$$
CREATE DEFINER=`dtuser`@`172.%` PROCEDURE `Select_LibraryFull` ()   BEGIN
	SELECT l.*, ss.secondary, st.experiment
	FROM Library l
	LEFT JOIN Structure_Secondary ss ON ss.LID=l.ID
	LEFT JOIN (SELECT sl.LID,s.SID, s.experiment FROM Structure s INNER JOIN Structure_Library sl ON sl.SID =s.SID GROUP BY sl.LID) st ON st.LID=l.ID
	ORDER BY contig;
END$$

DROP PROCEDURE IF EXISTS `Select_LofdControl`$$
CREATE DEFINER=`dtuser`@`172.%` PROCEDURE `Select_LofdControl` (IN `lids` TEXT, IN `types` VARCHAR(100))   BEGIN
	SELECT g.LID, g.read_index, g.`position`, g.predict_dwell as 'predict', sc.predict as 'control'
	FROM Lof g
	INNER JOIN Structure_Library sl ON sl.LID=g.LID 
	INNER JOIN Structure_Control sc ON sc.SID=sl.SID AND sc.`position`=g.`position`
	WHERE g.LID IN (lids) AND FIND_IN_SET(sc.type,`types`);
END$$

DROP PROCEDURE IF EXISTS `Select_LofdStructure`$$
CREATE DEFINER=`dtuser`@`172.%` PROCEDURE `Select_LofdStructure` (IN `lids` TEXT)   BEGIN
	SELECT g.LID, g.read_index, g.`position`, g.predict_dwell as 'predict', s.base_type, s.metric
	FROM Lof g
	INNER JOIN Structure_Library sl ON sl.LID=g.LID 
	INNER JOIN Structure s ON s.SID=sl.SID AND s.`position`=g.`position`
	WHERE g.LID IN (lids);
END$$

DROP PROCEDURE IF EXISTS `Select_LofsControl`$$
CREATE DEFINER=`dtuser`@`172.%` PROCEDURE `Select_LofsControl` (IN `lids` TEXT, IN `types` VARCHAR(100))   BEGIN
	SELECT g.LID, g.read_index, g.`position`, g.predict_signal as 'predict', sc.predict as 'control'
	FROM Lof g
	INNER JOIN Structure_Library sl ON sl.LID=g.LID 
	INNER JOIN Structure_Control sc ON sc.SID=sl.SID AND sc.`position`=g.`position`
	WHERE g.LID IN (lids) AND FIND_IN_SET(sc.type,`types`);
END$$

DROP PROCEDURE IF EXISTS `Select_LofsStructure`$$
CREATE DEFINER=`dtuser`@`172.%` PROCEDURE `Select_LofsStructure` (IN `lids` TEXT)   BEGIN
	SELECT g.LID, g.read_index, g.`position`, g.predict_signal as 'predict', s.base_type, s.metric
	FROM Lof g
	INNER JOIN Structure_Library sl ON sl.LID=g.LID 
	INNER JOIN Structure s ON s.SID=sl.SID AND s.`position`=g.`position`
	WHERE g.LID IN (lids);
END$$

DROP PROCEDURE IF EXISTS `Select_MaxClusters`$$
CREATE DEFINER=`dtuser`@`172.%` PROCEDURE `Select_MaxClusters` (IN `lids` TEXT)   BEGIN
	SELECT LID, c.contig, cluster, method, COUNT(*)/l.sequence_len as 'cluster_size' 
	FROM Clusters c 
	INNER JOIN Library l ON l.ID=c.LID
	WHERE LID IN (lids)
	GROUP BY cluster, method
	ORDER BY method, cluster;
END$$

DROP PROCEDURE IF EXISTS `Select_ModLID`$$
CREATE DEFINER=`dtuser`@`172.%` PROCEDURE `Select_ModLID` (IN `lids` TEXT)   BEGIN
	SELECT d.LID, d.contig, read_index, `position`, event_level_mean, event_length, l.run, l.temp, l.complex, l.type1, l.type2
	FROM Library l
	INNER JOIN Modified d ON d.LID=l.ID
	WHERE d.LID IN (`lids`) 
	ORDER BY contig, read_index, `position`, run, temp;
END$$

DROP PROCEDURE IF EXISTS `Select_PeaksUnmod`$$
CREATE DEFINER=`dtuser`@`172.%` PROCEDURE `Select_PeaksUnmod` (IN `mod_lids` TEXT, IN `unmod_lids` TEXT)   BEGIN
	DROP TEMPORARY TABLE IF EXISTS t1;

	CREATE TEMPORARY TABLE t1 AS
		SELECT contig, `position`, AVG(event_length) AS 'event_length_dmso', AVG(event_level_mean) AS 'event_level_mean_dmso'
		FROM Unmodified d
		WHERE d.LID IN (unmod_lids)
		GROUP BY contig, `position`;
	
	SELECT a.LID, a.contig, a.position, a.read_index, (ABS(event_length) - ABS(event_length_dmso)) as 'delta_dwell', (ABS(event_level_mean) - ABS(event_level_mean_dmso)) as 'delta_signal'   
	FROM Unmodified a
	INNER JOIN t1 ON t1.contig=a.contig AND t1.position=a.position
	WHERE a.LID IN (mod_lids)
	ORDER BY a.LID, a.contig, a.read_index, a.position;
END$$

DROP PROCEDURE IF EXISTS `Select_PutativeStructures`$$
CREATE DEFINER=`dtuser`@`172.%` PROCEDURE `Select_PutativeStructures` (IN `lids` TEXT)   BEGIN
	SELECT l.ID, l.contig, l.sequence, l.sequence_len, cs.secondary, cs.cluster, cs.mfe, cs.`method`
	FROM Library l 
	INNER JOIN Centroid_Secondary cs ON cs.LID=l.ID
	WHERE cs.LID IN (lids) AND cs.type='MFE' AND cs.mfe <> 0;
END$$

DROP PROCEDURE IF EXISTS `Select_RDAverage`$$
CREATE DEFINER=`dtuser`@`172.%` PROCEDURE `Select_RDAverage` (IN `lids` TEXT)   BEGIN
	SELECT * FROM Read_depth WHERE Read_depth.LID IN (lids);
END$$

DROP PROCEDURE IF EXISTS `Select_RdfContinue`$$
CREATE DEFINER=`dtuser`@`172.%` PROCEDURE `Select_RdfContinue` (IN `lids` TEXT, IN `lids2` TEXT)   BEGIN
	SELECT ssi.LID, ssi.read_index, 1 as 'completed' 
	FROM Structure_Secondary_Interactions ssi 
	INNER JOIN Structure_BPP sb ON sb.SSID =ssi.SSID
	INNER JOIN Structure_Probabilities sp ON sp.SSID =sb.SSID AND sp.BPPID =sb.BPPID
	WHERE ssi.LID IN (lids) AND ssi.LID2 IN (lids2)
	GROUP BY ssi.LID, ssi.read_index;
END$$

DROP PROCEDURE IF EXISTS `Select_RDFControl`$$
CREATE DEFINER=`dtuser`@`172.%` PROCEDURE `Select_RDFControl` (IN `lids` TEXT, IN `types` VARCHAR(100))   BEGIN
	SELECT g.LID, g.read_index, g.`position`, g.Predict as 'predict', sc.predict as 'control'
	FROM Read_depth_full g
	INNER JOIN Structure_Library sl ON sl.LID=g.LID 
	INNER JOIN Structure_Control sc ON sc.SID=sl.SID AND sc.`position`=g.`position`
	WHERE g.LID IN (lids) AND FIND_IN_SET(sc.type,`types`);
END$$

DROP PROCEDURE IF EXISTS `Select_RDFLID`$$
CREATE DEFINER=`dtuser`@`172.%` PROCEDURE `Select_RDFLID` (IN `lids` TEXT)   BEGIN
	SELECT * FROM Read_depth_full
	WHERE LID IN (lids);
END$$

DROP PROCEDURE IF EXISTS `Select_RDFStructure`$$
CREATE DEFINER=`dtuser`@`172.%` PROCEDURE `Select_RDFStructure` (IN `lids` TEXT)   BEGIN
	SELECT g.LID, g.read_index, g.`position`, g.Predict as 'predict', s.base_type, s.metric
	FROM Read_depth_full g
	INNER JOIN Structure_Library sl ON sl.LID=g.LID 
	INNER JOIN Structure s ON s.SID=sl.SID AND s.`position`=g.`position`
	WHERE g.LID IN (lids);
END$$

DROP PROCEDURE IF EXISTS `Select_ReadDepth`$$
CREATE DEFINER=`dtuser`@`172.%` PROCEDURE `Select_ReadDepth` (IN `lids` TEXT)   BEGIN
	SELECT LID, position, COUNT(read_index) as 'read_depth'
	FROM Unmodified u 
	WHERE LID IN (lids)
	GROUP BY LID, position
	UNION 
	SELECT LID, position, COUNT(read_index) as 'read_depth'
	FROM Modified m 
	WHERE LID IN (lids)
	GROUP BY LID, position;
END$$

DROP PROCEDURE IF EXISTS `Select_ReadDepthControl`$$
CREATE DEFINER=`dtuser`@`172.%` PROCEDURE `Select_ReadDepthControl` (IN `lids` TEXT, IN `types` VARCHAR(100))   BEGIN
	SELECT g.LID, g.`position`, g.Predict as 'predict', sc.predict as 'control'
	FROM Read_depth g
	INNER JOIN Structure_Library sl ON sl.LID=g.LID 
	INNER JOIN Structure_Control sc ON sc.SID=sl.SID AND sc.`position`=g.`position`
	WHERE g.LID IN (lids) AND FIND_IN_SET(sc.type,`types`);
END$$

DROP PROCEDURE IF EXISTS `Select_ReadDepthStructure`$$
CREATE DEFINER=`dtuser`@`172.%` PROCEDURE `Select_ReadDepthStructure` (IN `lids` TEXT)   BEGIN
	SELECT g.LID, g.`position`, g.Predict as 'predict', s.base_type, s.metric
	FROM Read_depth g
	INNER JOIN Structure_Library sl ON sl.LID=g.LID 
	INNER JOIN Structure s ON s.SID=sl.SID AND s.`position`=g.`position`
	WHERE g.LID IN (lids);
END$$

DROP PROCEDURE IF EXISTS `Select_SecondaryStructure`$$
CREATE DEFINER=`dtuser`@`172.%` PROCEDURE `Select_SecondaryStructure` (IN `lids` TEXT)   BEGIN
	SELECT s.ID, s.LID, s.contig, s.secondary, l.sequence, s.mfe
	FROM Structure_Secondary s
	INNER JOIN Library l ON l.ID=s.LID
	WHERE s.LID IN (lids);
END$$

DROP PROCEDURE IF EXISTS `Select_SecondaryStructures`$$
CREATE DEFINER=`dtuser`@`172.%` PROCEDURE `Select_SecondaryStructures` ()   BEGIN
	SELECT s.ID, s.LID, s.contig, s.secondary, l.sequence, mfe
	FROM Structure_Secondary s
	INNER JOIN Library l ON l.ID=s.LID;
END$$

DROP PROCEDURE IF EXISTS `Select_ShapeControl`$$
CREATE DEFINER=`dtuser`@`172.%` PROCEDURE `Select_ShapeControl` (IN `lids` TEXT, IN `type1` VARCHAR(100), IN `types` VARCHAR(100))   BEGIN
	SELECT sl.LID, g.`position`, g.predict, sc.predict as 'control', s.base_type, s.metric
	FROM Structure_Control g
	INNER JOIN Structure_Control sc ON sc.SID=g.SID AND sc.position=g.position
	INNER JOIN Structure_Library sl ON sl.SID=g.SID 
	INNER JOIN Structure s ON s.SID=sl.SID AND s.`position`=g.`position`
	WHERE sl.LID IN (lids) AND g.type=type1 AND FIND_IN_SET(sc.type,`types`);
END$$

DROP PROCEDURE IF EXISTS `Select_SignalControl`$$
CREATE DEFINER=`dtuser`@`172.%` PROCEDURE `Select_SignalControl` (IN `lids` TEXT, IN `types` VARCHAR(100))   BEGIN
	SELECT g.LID, g.read_index, g.`position`, g.predict_signal as 'predict', sc.predict as 'control'
	FROM Peaks g
	INNER JOIN Structure_Library sl ON sl.LID=g.LID 
	INNER JOIN Structure_Control sc ON sc.SID=sl.SID AND sc.`position`=g.`position`
	WHERE g.LID IN (lids) AND FIND_IN_SET(sc.type,`types`);
END$$

DROP PROCEDURE IF EXISTS `Select_SignalStructure`$$
CREATE DEFINER=`dtuser`@`172.%` PROCEDURE `Select_SignalStructure` (IN `lids` TEXT)   BEGIN
	SELECT g.LID, g.read_index, g.`position`, g.predict_signal as 'predict', s.base_type, s.metric
	FROM Peaks g
	INNER JOIN Structure_Library sl ON sl.LID=g.LID 
	INNER JOIN Structure s ON s.SID=sl.SID AND s.`position`=g.`position`
	WHERE g.LID IN (lids);
END$$

DROP PROCEDURE IF EXISTS `Select_UnmodLID`$$
CREATE DEFINER=`dtuser`@`172.%` PROCEDURE `Select_UnmodLID` (IN `lids` TEXT)   BEGIN
	SELECT d.LID, d.contig, read_index, `position`, event_level_mean, event_length, l.run, l.temp, l.complex, l.type1, l.type2
	FROM Library l
	INNER JOIN Unmodified d ON d.LID=l.ID
	WHERE d.LID IN (`lids`) 
	ORDER BY contig, read_index, `position`, run, temp;
END$$

DROP PROCEDURE IF EXISTS `Select_UnmodSSI`$$
CREATE DEFINER=`dtuser`@`172.%` PROCEDURE `Select_UnmodSSI` ()   BEGIN
	SELECT l.ID as 'LID'
	FROM Library l 
	WHERE l.ID NOT IN (SELECT l.ID
	FROM Library l 
	LEFT JOIN Structure_Secondary_Interactions ssi ON ssi.LID =l.ID
	WHERE l.is_modified = 0 AND ssi.read_index=-2
	GROUP BY l.ID) AND l.is_modified =0;
END$$

DROP PROCEDURE IF EXISTS `Select_ViennaControl`$$
CREATE DEFINER=`dtuser`@`172.%` PROCEDURE `Select_ViennaControl` (IN `lids` TEXT, IN `types` VARCHAR(100))   BEGIN
	SELECT g.LID, g.read_index, g.`position`, g.Base_pair_prob as 'predict', sc.predict as 'control'
	FROM Read_depth_full g
	INNER JOIN Structure_Library sl ON sl.LID=g.LID 
	INNER JOIN Structure_Control sc ON sc.SID=sl.SID AND sc.`position`=g.`position`
	WHERE g.LID IN (lids) AND FIND_IN_SET(sc.type,`types`);
END$$

DROP PROCEDURE IF EXISTS `Select_ViennaStructure`$$
CREATE DEFINER=`dtuser`@`172.%` PROCEDURE `Select_ViennaStructure` (IN `lids` TEXT)   BEGIN
	SELECT g.LID, g.read_index, g.`position`, g.Base_pair_prob as 'predict', s.base_type, s.metric
	FROM Read_depth_full g
	INNER JOIN Structure_Library sl ON sl.LID=g.LID 
	INNER JOIN Structure s ON s.SID=sl.SID AND s.`position`=g.`position`
	WHERE g.LID IN (lids);
END$$

DROP PROCEDURE IF EXISTS `Update_BPP`$$
CREATE DEFINER=`dtuser`@`172.%` PROCEDURE `Update_BPP` (IN `lids` TEXT, IN `threshold` FLOAT, IN `rx_threshold` FLOAT)   BEGIN
	UPDATE Read_depth_full rdf 
	INNER JOIN 
	(SELECT LID, base, read_index, MAX(max_prob) AS max_prob
	FROM (
	    SELECT 
	        rdf.LID, 
	        sp.base1 - 1 AS base, 
	        ssi.read_index, 
	        MAX(sp.probability) AS max_prob
	    FROM Read_depth_full rdf
	    JOIN Structure_Secondary_Interactions ssi 
	        ON ssi.LID2 = rdf.LID AND ssi.read_index = rdf.read_index
	    JOIN Structure_BPP sb 
	        ON sb.SSID = ssi.SSID
	    JOIN Structure_Probabilities sp 
	        ON sp.SSID = ssi.SSID
	    WHERE rdf.LID IN (`lids`) 
	        AND ssi.type = 'MFE' 
	        AND sb.type_interaction IN ('A', 'AA')
	    GROUP BY rdf.LID, sp.base1, ssi.read_index
	
	    UNION
	
	    SELECT 
	        rdf.LID, 
	        sp.base2 - 1 AS base, 
	        ssi.read_index, 
	        MAX(sp.probability) AS max_prob
	    FROM Read_depth_full rdf
	    JOIN Structure_Secondary_Interactions ssi 
	        ON ssi.LID2 = rdf.LID AND ssi.read_index = rdf.read_index
	    JOIN Structure_BPP sb 
	        ON sb.SSID = ssi.SSID
	    JOIN Structure_Probabilities sp 
	        ON sp.SSID = ssi.SSID
	    WHERE rdf.LID IN (`lids`) 
	        AND ssi.type = 'MFE' 
	        AND sb.type_interaction IN ('A', 'AA')
	    GROUP BY rdf.LID, sp.base2, ssi.read_index
	) AS C
	GROUP BY base, read_index) AS A ON rdf.LID=A.LID AND rdf.read_index=A.read_index AND rdf.position = A.base
	SET rdf.base_pair_prob = A.max_prob, rdf.completed=1
	WHERE rdf.LID IN (`lids`);
	    
	UPDATE Read_depth_full rdf 
	SET rdf.Predict=1
	WHERE rdf.LID IN (`lids`)  AND rdf.base_pair_prob > `threshold` AND rdf.Reactivity_score < `rx_threshold`;

END$$

DELIMITER ;

-- --------------------------------------------------------

--
-- Table structure for table `Basecall`
--

DROP TABLE IF EXISTS `Basecall`;
CREATE TABLE `Basecall` (
  `LID` int(11) NOT NULL DEFAULT 0,
  `position` int(11) NOT NULL,
  `contig` varchar(255) NOT NULL,
  `insertion` float NOT NULL DEFAULT 0,
  `mismatch` float NOT NULL DEFAULT 0,
  `deletion` float NOT NULL DEFAULT 0,
  `quality` float NOT NULL DEFAULT 0,
  `basecall_reactivity` float NOT NULL DEFAULT 0,
  `aligned_reads` int(11) NOT NULL DEFAULT 0,
  `sequence` char(1) DEFAULT NULL,
  `timestamp` datetime NOT NULL DEFAULT current_timestamp(),
  `ID` int(10) UNSIGNED NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_general_ci;

-- --------------------------------------------------------

--
-- Table structure for table `Basecall_Peaks`
--

DROP TABLE IF EXISTS `Basecall_Peaks`;
CREATE TABLE `Basecall_Peaks` (
  `LID` int(11) NOT NULL DEFAULT 0,
  `position` int(11) NOT NULL,
  `contig` varchar(255) NOT NULL,
  `is_peak` int(11) NOT NULL,
  `peak_height` float DEFAULT NULL,
  `insertion` float DEFAULT NULL,
  `mismatch` float DEFAULT NULL,
  `deletion` float NOT NULL DEFAULT 0,
  `quality` float NOT NULL DEFAULT 0,
  `basecall_reactivity` float NOT NULL DEFAULT 0,
  `aligned_reads` int(11) NOT NULL DEFAULT 0,
  `predict` int(11) NOT NULL,
  `timestamp` datetime NOT NULL,
  `varna` int(11) DEFAULT NULL,
  `ID` int(10) UNSIGNED NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_general_ci;

-- --------------------------------------------------------

--
-- Table structure for table `Centroids`
--

DROP TABLE IF EXISTS `Centroids`;
CREATE TABLE `Centroids` (
  `ID` int(11) NOT NULL,
  `LID` int(11) NOT NULL DEFAULT 0,
  `contig` varchar(255) NOT NULL,
  `position` int(11) NOT NULL,
  `method` varchar(100) NOT NULL,
  `cluster` int(11) NOT NULL,
  `centroid` float NOT NULL,
  `timestamp` datetime NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb3 COLLATE=utf8mb3_unicode_ci;

-- --------------------------------------------------------

--
-- Table structure for table `Centroid_BPP`
--

DROP TABLE IF EXISTS `Centroid_BPP`;
CREATE TABLE `Centroid_BPP` (
  `ID` int(11) NOT NULL,
  `SSID` bigint(20) NOT NULL,
  `BPPID` bigint(20) NOT NULL,
  `type_interaction` varchar(100) NOT NULL,
  `mfe` float DEFAULT NULL,
  `Timestamp` datetime NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb3 COLLATE=utf8mb3_unicode_ci;

-- --------------------------------------------------------

--
-- Table structure for table `Centroid_Distance`
--

DROP TABLE IF EXISTS `Centroid_Distance`;
CREATE TABLE `Centroid_Distance` (
  `ID` int(11) NOT NULL,
  `LID` int(11) NOT NULL,
  `cluster` int(11) NOT NULL,
  `distance` float DEFAULT NULL,
  `method` varchar(100) NOT NULL,
  `timestamp` datetime NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb3 COLLATE=utf8mb3_unicode_ci;

-- --------------------------------------------------------

--
-- Table structure for table `Centroid_Probabilities`
--

DROP TABLE IF EXISTS `Centroid_Probabilities`;
CREATE TABLE `Centroid_Probabilities` (
  `ID` int(11) NOT NULL,
  `SSID` bigint(20) NOT NULL,
  `BPPID` bigint(20) NOT NULL,
  `base1` int(11) NOT NULL,
  `base2` int(11) NOT NULL,
  `probability` float NOT NULL,
  `Timestamp` datetime NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb3 COLLATE=utf8mb3_unicode_ci;

-- --------------------------------------------------------

--
-- Table structure for table `Centroid_Regions`
--

DROP TABLE IF EXISTS `Centroid_Regions`;
CREATE TABLE `Centroid_Regions` (
  `ID` int(11) NOT NULL,
  `LID` int(11) NOT NULL,
  `contig` varchar(100) NOT NULL,
  `position` int(11) NOT NULL,
  `is_conserved` int(11) NOT NULL DEFAULT 0,
  `method` varchar(100) NOT NULL,
  `timestamp` datetime NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb3 COLLATE=utf8mb3_unicode_ci;

-- --------------------------------------------------------

--
-- Table structure for table `Centroid_Secondary`
--

DROP TABLE IF EXISTS `Centroid_Secondary`;
CREATE TABLE `Centroid_Secondary` (
  `ID` int(11) NOT NULL,
  `LID` int(11) NOT NULL,
  `cluster` int(11) NOT NULL,
  `contig` varchar(100) NOT NULL,
  `secondary` text NOT NULL,
  `type` varchar(100) DEFAULT NULL,
  `mfe` float DEFAULT NULL,
  `frequency` float DEFAULT NULL,
  `diversity` float DEFAULT NULL,
  `distance` float DEFAULT NULL,
  `mea` float DEFAULT NULL,
  `method` varchar(100) DEFAULT NULL,
  `Timestamp` datetime NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb3 COLLATE=utf8mb3_unicode_ci;

-- --------------------------------------------------------

--
-- Table structure for table `Centroid_Secondary_Interactions`
--

DROP TABLE IF EXISTS `Centroid_Secondary_Interactions`;
CREATE TABLE `Centroid_Secondary_Interactions` (
  `ID` int(11) NOT NULL,
  `SSID` bigint(20) NOT NULL DEFAULT 0,
  `LID1` int(11) NOT NULL,
  `cluster1` int(11) NOT NULL,
  `contig1` varchar(100) NOT NULL,
  `LID2` int(11) NOT NULL,
  `cluster2` int(11) NOT NULL,
  `contig2` varchar(100) NOT NULL,
  `secondary` text NOT NULL,
  `mfe` float DEFAULT NULL,
  `frequency` float DEFAULT NULL,
  `diversity` float DEFAULT NULL,
  `distance` float DEFAULT NULL,
  `mea` float DEFAULT NULL,
  `deltag` float DEFAULT NULL,
  `type` varchar(100) DEFAULT NULL,
  `method` varchar(100) DEFAULT NULL,
  `Timestamp` datetime NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb3 COLLATE=utf8mb3_unicode_ci;

-- --------------------------------------------------------

--
-- Table structure for table `Clusters`
--

DROP TABLE IF EXISTS `Clusters`;
CREATE TABLE `Clusters` (
  `ID` int(11) NOT NULL,
  `LID` int(11) NOT NULL DEFAULT 0,
  `contig` varchar(255) NOT NULL,
  `read_index` int(11) NOT NULL,
  `cluster` int(11) NOT NULL,
  `position` int(11) NOT NULL,
  `sequence` char(1) DEFAULT NULL,
  `Reactivity` float NOT NULL,
  `Predict` int(11) DEFAULT NULL,
  `timestamp` datetime NOT NULL,
  `method` varchar(100) DEFAULT 'kmeans'
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb3 COLLATE=utf8mb3_unicode_ci;

-- --------------------------------------------------------

--
-- Table structure for table `Gmm`
--

DROP TABLE IF EXISTS `Gmm`;
CREATE TABLE `Gmm` (
  `LID` int(11) NOT NULL DEFAULT 0,
  `position` int(11) NOT NULL,
  `contig` varchar(255) NOT NULL,
  `read_depth` int(11) DEFAULT NULL,
  `percent_modified` float DEFAULT NULL,
  `predict` int(11) NOT NULL,
  `timestamp` datetime DEFAULT NULL,
  `ID` int(10) UNSIGNED NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_general_ci;

-- --------------------------------------------------------

--
-- Table structure for table `Library`
--

DROP TABLE IF EXISTS `Library`;
CREATE TABLE `Library` (
  `ID` int(11) NOT NULL,
  `contig` varchar(255) NOT NULL DEFAULT 'c',
  `sequence_name` varchar(255) NOT NULL,
  `sequence_len` int(11) NOT NULL,
  `temp` float NOT NULL DEFAULT 37,
  `is_modified` tinyint(1) NOT NULL DEFAULT 1,
  `type1` varchar(100) NOT NULL DEFAULT 'dmso',
  `type2` varchar(100) NOT NULL DEFAULT 'acim',
  `complex` int(11) NOT NULL DEFAULT 0,
  `run` int(11) NOT NULL DEFAULT 1,
  `sequence` text NOT NULL,
  `timestamp` datetime NOT NULL DEFAULT current_timestamp(),
  `is_putative` binary(1) NOT NULL DEFAULT '0'
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb3 COLLATE=utf8mb3_unicode_ci;

-- --------------------------------------------------------

--
-- Table structure for table `Lof`
--

DROP TABLE IF EXISTS `Lof`;
CREATE TABLE `Lof` (
  `LID` int(11) NOT NULL DEFAULT 0,
  `position` int(11) NOT NULL,
  `contig` varchar(255) NOT NULL,
  `read_index` int(11) NOT NULL,
  `predict_dwell` int(11) DEFAULT NULL,
  `varna_dwell` int(11) DEFAULT NULL,
  `predict_signal` int(11) DEFAULT NULL,
  `varna_signal` int(11) DEFAULT NULL,
  `timestamp` datetime DEFAULT NULL,
  `ID` int(10) UNSIGNED NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_general_ci;

-- --------------------------------------------------------

--
-- Table structure for table `Modified`
--

DROP TABLE IF EXISTS `Modified`;
CREATE TABLE `Modified` (
  `LID` int(11) NOT NULL DEFAULT 0,
  `ID` int(11) NOT NULL,
  `contig` varchar(255) NOT NULL,
  `position` int(11) NOT NULL,
  `reference_kmer` varchar(10) NOT NULL,
  `read_index` int(11) NOT NULL,
  `event_level_mean` float NOT NULL,
  `event_length` float NOT NULL,
  `event_stdv` float NOT NULL,
  `event_level_mean_norm` float DEFAULT NULL,
  `event_length_norm` float DEFAULT NULL,
  `event_level_mean_raw` float DEFAULT NULL,
  `event_length_raw` float DEFAULT NULL,
  `timestamp` datetime NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb3 COLLATE=utf8mb3_unicode_ci;

-- --------------------------------------------------------

--
-- Table structure for table `Peaks`
--

DROP TABLE IF EXISTS `Peaks`;
CREATE TABLE `Peaks` (
  `LID` int(11) NOT NULL DEFAULT 0,
  `read_index` int(11) NOT NULL,
  `position` int(11) NOT NULL,
  `contig` varchar(255) NOT NULL,
  `delta_signal` float NOT NULL,
  `predict_signal` int(11) NOT NULL,
  `varna_signal` int(11) DEFAULT NULL,
  `delta_dwell` float NOT NULL,
  `predict_dwell` int(11) NOT NULL,
  `varna_dwell` int(11) DEFAULT NULL,
  `timestamp` datetime DEFAULT NULL,
  `ID` int(10) UNSIGNED NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_general_ci;

-- --------------------------------------------------------

--
-- Table structure for table `RAW`
--

DROP TABLE IF EXISTS `RAW`;
CREATE TABLE `RAW` (
  `LID` int(11) NOT NULL DEFAULT 0,
  `ID` int(11) NOT NULL,
  `contig` varchar(255) NOT NULL,
  `position` int(11) NOT NULL,
  `reference_kmer` varchar(10) NOT NULL,
  `read_index` int(11) NOT NULL,
  `event_level_mean` float NOT NULL,
  `event_length` float NOT NULL,
  `event_stdv` float NOT NULL,
  `event_level_mean_norm` float DEFAULT NULL,
  `event_length_norm` float DEFAULT NULL,
  `event_level_mean_std` float DEFAULT NULL,
  `event_length_std` float DEFAULT NULL,
  `temp` int(11) DEFAULT 37,
  `run` int(11) DEFAULT 1,
  `timestamp` datetime NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb3 COLLATE=utf8mb3_unicode_ci;

-- --------------------------------------------------------

--
-- Table structure for table `RAW_MOD`
--

DROP TABLE IF EXISTS `RAW_MOD`;
CREATE TABLE `RAW_MOD` (
  `LID` int(11) NOT NULL DEFAULT 0,
  `ID` int(11) NOT NULL,
  `contig` varchar(255) NOT NULL,
  `position` int(11) NOT NULL,
  `reference_kmer` varchar(10) NOT NULL,
  `read_index` int(11) NOT NULL,
  `event_level_mean` float NOT NULL,
  `event_length` float NOT NULL,
  `event_stdv` float NOT NULL,
  `event_level_mean_norm` float DEFAULT NULL,
  `event_length_norm` float DEFAULT NULL,
  `event_level_mean_std` float DEFAULT NULL,
  `event_length_std` float DEFAULT NULL,
  `temp` int(11) DEFAULT 37,
  `run` int(11) DEFAULT 1,
  `timestamp` datetime NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb3 COLLATE=utf8mb3_unicode_ci;

-- --------------------------------------------------------

--
-- Table structure for table `Reactivity_full`
--

DROP TABLE IF EXISTS `Reactivity_full`;
CREATE TABLE `Reactivity_full` (
  `ID` int(11) NOT NULL,
  `LID` int(11) NOT NULL DEFAULT 0,
  `contig` varchar(255) NOT NULL,
  `read_index` int(11) NOT NULL,
  `position` int(11) NOT NULL,
  `Sequence` char(1) NOT NULL,
  `Predict_BC` int(11) NOT NULL,
  `Predict_Signal` int(11) NOT NULL,
  `Predict_Dwell` int(11) NOT NULL,
  `Predict_Lofs` int(11) NOT NULL,
  `Predict_Lofd` int(11) NOT NULL,
  `Predict_Gmm` int(11) NOT NULL,
  `Reactivity` int(11) NOT NULL,
  `Predict` int(11) NOT NULL,
  `Varna` int(11) NOT NULL,
  `timestamp` datetime NOT NULL,
  `bae_pair_prob` float DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb3 COLLATE=utf8mb3_unicode_ci;

-- --------------------------------------------------------

--
-- Table structure for table `Read_depth`
--

DROP TABLE IF EXISTS `Read_depth`;
CREATE TABLE `Read_depth` (
  `ID` int(11) NOT NULL,
  `LID` int(11) NOT NULL DEFAULT 0,
  `contig` varchar(255) NOT NULL,
  `read_depth` int(11) NOT NULL,
  `position` int(11) NOT NULL,
  `read_index` int(11) NOT NULL DEFAULT -2,
  `Out_num` int(11) NOT NULL,
  `In_num` int(11) NOT NULL,
  `Percent_modified` float NOT NULL,
  `Rnafold_shape_reactivity` float NOT NULL,
  `Base_pair_prob` float NOT NULL DEFAULT 0,
  `Predict` int(11) NOT NULL,
  `Varna` int(11) NOT NULL,
  `timestamp` datetime NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb3 COLLATE=utf8mb3_unicode_ci;

-- --------------------------------------------------------

--
-- Table structure for table `Read_depth_full`
--

DROP TABLE IF EXISTS `Read_depth_full`;
CREATE TABLE `Read_depth_full` (
  `ID` int(11) NOT NULL,
  `LID` int(11) NOT NULL DEFAULT 0,
  `contig` varchar(255) NOT NULL,
  `read_index` int(11) NOT NULL,
  `position` int(11) NOT NULL,
  `Sequence` char(1) NOT NULL,
  `Predict_BC` int(11) NOT NULL,
  `Predict_Signal` int(11) NOT NULL,
  `Predict_Dwell` int(11) NOT NULL,
  `Predict_Lofs` int(11) NOT NULL,
  `Predict_Lofd` int(11) NOT NULL,
  `Predict_Gmm` int(11) NOT NULL,
  `Reactivity` int(11) NOT NULL,
  `Reactivity_score` float NOT NULL,
  `Predict` int(11) NOT NULL,
  `Varna` int(11) NOT NULL,
  `timestamp` datetime NOT NULL,
  `base_pair_prob` float DEFAULT NULL,
  `completed` tinyint(1) NOT NULL DEFAULT 0
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb3 COLLATE=utf8mb3_unicode_ci;

-- --------------------------------------------------------

--
-- Table structure for table `Structure`
--

DROP TABLE IF EXISTS `Structure`;
CREATE TABLE `Structure` (
  `ID` int(10) UNSIGNED NOT NULL,
  `SID` int(11) NOT NULL DEFAULT 0,
  `position` int(11) NOT NULL,
  `contig` varchar(255) NOT NULL,
  `sequence` char(1) NOT NULL,
  `base_type` char(1) DEFAULT NULL,
  `structure_type` char(10) DEFAULT NULL,
  `experiment` varchar(255) NOT NULL DEFAULT 'xray crystallography',
  `timestamp` datetime DEFAULT current_timestamp(),
  `metric` varchar(100) DEFAULT NULL,
  `metric_num` varchar(100) DEFAULT NULL,
  `accessibility` varchar(100) NOT NULL DEFAULT 'Y'
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_general_ci;

-- --------------------------------------------------------

--
-- Table structure for table `Structure_BPP`
--

DROP TABLE IF EXISTS `Structure_BPP`;
CREATE TABLE `Structure_BPP` (
  `ID` int(11) NOT NULL,
  `SSID` bigint(20) NOT NULL,
  `BPPID` bigint(20) NOT NULL,
  `type_interaction` varchar(100) NOT NULL,
  `mfe` float DEFAULT NULL,
  `Timestamp` datetime NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb3 COLLATE=utf8mb3_unicode_ci;

-- --------------------------------------------------------

--
-- Table structure for table `Structure_Control`
--

DROP TABLE IF EXISTS `Structure_Control`;
CREATE TABLE `Structure_Control` (
  `SID` int(11) NOT NULL DEFAULT 0,
  `position` int(11) NOT NULL,
  `contig` varchar(255) NOT NULL,
  `type` varchar(10) DEFAULT NULL,
  `reactivity` float NOT NULL,
  `sequence` char(1) NOT NULL,
  `predict` int(11) NOT NULL,
  `timestamp` timestamp NOT NULL DEFAULT current_timestamp(),
  `ID` int(10) UNSIGNED NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_general_ci;

-- --------------------------------------------------------

--
-- Table structure for table `Structure_Library`
--

DROP TABLE IF EXISTS `Structure_Library`;
CREATE TABLE `Structure_Library` (
  `ID` int(11) NOT NULL,
  `LID` int(11) NOT NULL,
  `SID` int(11) NOT NULL DEFAULT 0,
  `contig` varchar(100) DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb3 COLLATE=utf8mb3_unicode_ci;

-- --------------------------------------------------------

--
-- Table structure for table `Structure_Metric`
--

DROP TABLE IF EXISTS `Structure_Metric`;
CREATE TABLE `Structure_Metric` (
  `ID` int(11) NOT NULL,
  `LID` int(11) NOT NULL DEFAULT 0,
  `read_index` int(11) NOT NULL,
  `algorithm` varchar(100) DEFAULT NULL,
  `metric` varchar(100) NOT NULL,
  `threshold` float NOT NULL DEFAULT 0,
  `ppv` float DEFAULT NULL,
  `accuracy` float DEFAULT NULL,
  `sensitivity` float NOT NULL DEFAULT 0,
  `specificity` float NOT NULL DEFAULT 0,
  `tp` int(11) DEFAULT NULL,
  `fp` int(11) DEFAULT NULL,
  `tn` int(11) DEFAULT NULL,
  `fn` int(11) DEFAULT NULL,
  `pearsons` float DEFAULT NULL,
  `mannwhit` float NOT NULL,
  `aggreement` float NOT NULL,
  `timestamp` datetime DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb3 COLLATE=utf8mb3_unicode_ci;

-- --------------------------------------------------------

--
-- Table structure for table `Structure_Metric_Secondary`
--

DROP TABLE IF EXISTS `Structure_Metric_Secondary`;
CREATE TABLE `Structure_Metric_Secondary` (
  `ID` int(11) NOT NULL,
  `LID` int(11) NOT NULL DEFAULT 0,
  `read_index` int(11) NOT NULL,
  `algorithm` varchar(100) DEFAULT NULL,
  `metric` varchar(100) NOT NULL,
  `number_total` int(11) NOT NULL DEFAULT 0,
  `number_detected` int(11) NOT NULL DEFAULT 0,
  `bases_covered` int(11) NOT NULL DEFAULT 0,
  `bases_total` int(11) NOT NULL DEFAULT 0,
  `threshold` float NOT NULL DEFAULT 0,
  `timestamp` datetime DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb3 COLLATE=utf8mb3_unicode_ci;

-- --------------------------------------------------------

--
-- Table structure for table `Structure_Probabilities`
--

DROP TABLE IF EXISTS `Structure_Probabilities`;
CREATE TABLE `Structure_Probabilities` (
  `ID` int(11) NOT NULL,
  `SSID` bigint(20) NOT NULL,
  `BPPID` bigint(20) NOT NULL,
  `base1` int(11) NOT NULL,
  `base2` int(11) NOT NULL,
  `probability` float NOT NULL,
  `Timestamp` datetime NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb3 COLLATE=utf8mb3_unicode_ci;

-- --------------------------------------------------------

--
-- Table structure for table `Structure_Secondary`
--

DROP TABLE IF EXISTS `Structure_Secondary`;
CREATE TABLE `Structure_Secondary` (
  `ID` int(11) NOT NULL,
  `LID` int(11) NOT NULL,
  `contig` varchar(100) NOT NULL,
  `secondary` text NOT NULL,
  `timestamp` datetime NOT NULL,
  `mfe` float DEFAULT NULL,
  `psuedoknot` text DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb3 COLLATE=utf8mb3_unicode_ci;

-- --------------------------------------------------------

--
-- Table structure for table `Structure_Secondary_Interactions`
--

DROP TABLE IF EXISTS `Structure_Secondary_Interactions`;
CREATE TABLE `Structure_Secondary_Interactions` (
  `ID` int(11) NOT NULL,
  `SSID` bigint(20) NOT NULL,
  `LID` int(11) NOT NULL,
  `contig` varchar(100) NOT NULL,
  `LID2` int(11) NOT NULL,
  `contig2` varchar(100) NOT NULL,
  `read_index` int(11) NOT NULL DEFAULT -2,
  `secondary` text NOT NULL,
  `type` varchar(100) NOT NULL,
  `mfe` float DEFAULT NULL,
  `frequency` float DEFAULT NULL,
  `deltag` float DEFAULT NULL,
  `mea` float DEFAULT NULL,
  `Timestamp` datetime NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb3 COLLATE=utf8mb3_unicode_ci;

-- --------------------------------------------------------

--
-- Table structure for table `Unmodified`
--

DROP TABLE IF EXISTS `Unmodified`;
CREATE TABLE `Unmodified` (
  `LID` int(11) NOT NULL DEFAULT 0,
  `ID` int(11) NOT NULL,
  `contig` varchar(255) NOT NULL,
  `position` int(11) NOT NULL,
  `reference_kmer` varchar(10) NOT NULL,
  `read_index` int(11) NOT NULL,
  `event_level_mean` float NOT NULL,
  `event_length` float NOT NULL,
  `event_stdv` float NOT NULL,
  `event_level_mean_norm` float DEFAULT NULL,
  `event_length_norm` float DEFAULT NULL,
  `event_level_mean_raw` float DEFAULT NULL,
  `event_length_raw` float DEFAULT NULL,
  `timestamp` datetime NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb3 COLLATE=utf8mb3_unicode_ci;

--
-- Indexes for dumped tables
--

--
-- Indexes for table `Basecall`
--
ALTER TABLE `Basecall`
  ADD PRIMARY KEY (`ID`);

--
-- Indexes for table `Basecall_Peaks`
--
ALTER TABLE `Basecall_Peaks`
  ADD PRIMARY KEY (`ID`),
  ADD UNIQUE KEY `b1peak` (`LID`,`position`) USING BTREE;

--
-- Indexes for table `Centroids`
--
ALTER TABLE `Centroids`
  ADD PRIMARY KEY (`ID`),
  ADD UNIQUE KEY `Centroids_LID_IDX` (`LID`,`contig`,`position`,`cluster`,`method`) USING BTREE;

--
-- Indexes for table `Centroid_BPP`
--
ALTER TABLE `Centroid_BPP`
  ADD PRIMARY KEY (`ID`),
  ADD KEY `Structure_BPP_SSID_IDX` (`type_interaction`) USING BTREE,
  ADD KEY `maxtype` (`type_interaction`) USING BTREE,
  ADD KEY `Structure_BPP_SSID` (`SSID`,`BPPID`) USING BTREE;

--
-- Indexes for table `Centroid_Distance`
--
ALTER TABLE `Centroid_Distance`
  ADD PRIMARY KEY (`ID`),
  ADD UNIQUE KEY `Centroid_Distance_LID_IDX` (`LID`,`cluster`,`method`) USING BTREE;

--
-- Indexes for table `Centroid_Probabilities`
--
ALTER TABLE `Centroid_Probabilities`
  ADD PRIMARY KEY (`ID`),
  ADD KEY `max_prob1` (`base1`),
  ADD KEY `max_prob2` (`base2`),
  ADD KEY `max_prob3` (`base1`,`base2`),
  ADD KEY `maxbase1` (`base1`),
  ADD KEY `maxbase2` (`base2`),
  ADD KEY `Structure_Probabilities_SSID_IDX` (`SSID`,`BPPID`) USING BTREE,
  ADD KEY `Structure_Probabilities_Full` (`SSID`,`BPPID`,`base1`,`base2`,`probability`) USING BTREE;

--
-- Indexes for table `Centroid_Regions`
--
ALTER TABLE `Centroid_Regions`
  ADD PRIMARY KEY (`ID`);

--
-- Indexes for table `Centroid_Secondary`
--
ALTER TABLE `Centroid_Secondary`
  ADD PRIMARY KEY (`ID`);

--
-- Indexes for table `Centroid_Secondary_Interactions`
--
ALTER TABLE `Centroid_Secondary_Interactions`
  ADD PRIMARY KEY (`ID`);

--
-- Indexes for table `Clusters`
--
ALTER TABLE `Clusters`
  ADD PRIMARY KEY (`ID`),
  ADD UNIQUE KEY `Clusters_LID_IDX` (`LID`,`contig`,`read_index`,`position`,`cluster`,`method`) USING BTREE;

--
-- Indexes for table `Gmm`
--
ALTER TABLE `Gmm`
  ADD PRIMARY KEY (`ID`),
  ADD UNIQUE KEY `g1` (`LID`,`position`) USING BTREE;

--
-- Indexes for table `Library`
--
ALTER TABLE `Library`
  ADD PRIMARY KEY (`ID`),
  ADD UNIQUE KEY `sequence_name` (`sequence_name`,`temp`,`type1`,`type2`,`complex`,`run`) USING BTREE;

--
-- Indexes for table `Lof`
--
ALTER TABLE `Lof`
  ADD PRIMARY KEY (`ID`),
  ADD UNIQUE KEY `lofs` (`LID`,`read_index`,`position`,`ID`) USING BTREE,
  ADD KEY `lofs_control` (`LID`,`position`,`predict_signal`),
  ADD KEY `lofd_control` (`LID`,`position`,`predict_dwell`);

--
-- Indexes for table `Modified`
--
ALTER TABLE `Modified`
  ADD PRIMARY KEY (`ID`),
  ADD UNIQUE KEY `acim_idx_temp` (`LID`,`read_index`,`position`,`ID`) USING BTREE,
  ADD KEY `select_mod_idx` (`LID`,`contig`,`read_index`,`position`);

--
-- Indexes for table `Peaks`
--
ALTER TABLE `Peaks`
  ADD PRIMARY KEY (`ID`),
  ADD UNIQUE KEY `signal1` (`LID`,`read_index`,`position`,`ID`) USING BTREE,
  ADD KEY `signal_dwell` (`LID`,`position`,`predict_dwell`,`ID`) USING BTREE;

--
-- Indexes for table `RAW`
--
ALTER TABLE `RAW`
  ADD PRIMARY KEY (`ID`),
  ADD UNIQUE KEY `dmso_idx_temp` (`LID`,`read_index`,`position`) USING BTREE,
  ADD KEY `Dmso_contig_IDX` (`LID`,`position`) USING BTREE,
  ADD KEY `select_unmod_indx` (`LID`,`contig`,`read_index`,`position`);

--
-- Indexes for table `RAW_MOD`
--
ALTER TABLE `RAW_MOD`
  ADD PRIMARY KEY (`ID`),
  ADD UNIQUE KEY `dmso_idx_temp` (`LID`,`read_index`,`position`) USING BTREE,
  ADD KEY `Dmso_contig_IDX` (`LID`,`position`) USING BTREE,
  ADD KEY `select_unmod_indx` (`LID`,`contig`,`read_index`,`position`);

--
-- Indexes for table `Reactivity_full`
--
ALTER TABLE `Reactivity_full`
  ADD PRIMARY KEY (`ID`);

--
-- Indexes for table `Read_depth`
--
ALTER TABLE `Read_depth`
  ADD PRIMARY KEY (`ID`),
  ADD KEY `rd_control` (`LID`,`position`,`Predict`);

--
-- Indexes for table `Read_depth_full`
--
ALTER TABLE `Read_depth_full`
  ADD PRIMARY KEY (`ID`),
  ADD KEY `rdf_control` (`LID`,`position`,`Predict`),
  ADD KEY `rdf_continue` (`LID`,`read_index`),
  ADD KEY `rdf_prob` (`LID`,`read_index`,`position`) USING BTREE;

--
-- Indexes for table `Structure`
--
ALTER TABLE `Structure`
  ADD PRIMARY KEY (`ID`);

--
-- Indexes for table `Structure_BPP`
--
ALTER TABLE `Structure_BPP`
  ADD PRIMARY KEY (`ID`),
  ADD KEY `maxtype` (`type_interaction`) USING BTREE,
  ADD KEY `Structure_BPP_SSID` (`SSID`,`BPPID`) USING BTREE,
  ADD KEY `ssid_type` (`SSID`,`type_interaction`) USING BTREE;

--
-- Indexes for table `Structure_Control`
--
ALTER TABLE `Structure_Control`
  ADD PRIMARY KEY (`ID`);

--
-- Indexes for table `Structure_Library`
--
ALTER TABLE `Structure_Library`
  ADD PRIMARY KEY (`ID`);

--
-- Indexes for table `Structure_Metric`
--
ALTER TABLE `Structure_Metric`
  ADD PRIMARY KEY (`ID`);

--
-- Indexes for table `Structure_Metric_Secondary`
--
ALTER TABLE `Structure_Metric_Secondary`
  ADD PRIMARY KEY (`ID`);

--
-- Indexes for table `Structure_Probabilities`
--
ALTER TABLE `Structure_Probabilities`
  ADD PRIMARY KEY (`ID`),
  ADD KEY `max_prob1` (`base1`),
  ADD KEY `max_prob2` (`base2`),
  ADD KEY `max_prob3` (`base1`,`base2`),
  ADD KEY `maxbase1` (`base1`),
  ADD KEY `maxbase2` (`base2`),
  ADD KEY `Structure_Probabilities_SSID_IDX` (`SSID`,`BPPID`) USING BTREE,
  ADD KEY `Structure_Probabilities_Full` (`SSID`,`BPPID`,`base1`,`base2`,`probability`) USING BTREE,
  ADD KEY `sp_rdf` (`SSID`,`base1`) USING BTREE,
  ADD KEY `sp_rdf2` (`SSID`,`base2`) USING BTREE;

--
-- Indexes for table `Structure_Secondary`
--
ALTER TABLE `Structure_Secondary`
  ADD PRIMARY KEY (`ID`);

--
-- Indexes for table `Structure_Secondary_Interactions`
--
ALTER TABLE `Structure_Secondary_Interactions`
  ADD PRIMARY KEY (`ID`),
  ADD KEY `max_prob0` (`LID`,`read_index`),
  ADD KEY `max_prob00` (`read_index`),
  ADD KEY `lids` (`LID`,`LID2`) USING BTREE,
  ADD KEY `ssid` (`SSID`,`LID`,`LID2`) USING BTREE,
  ADD KEY `ssi_rdf` (`LID2`,`read_index`) USING BTREE;

--
-- Indexes for table `Unmodified`
--
ALTER TABLE `Unmodified`
  ADD PRIMARY KEY (`ID`),
  ADD UNIQUE KEY `dmso_idx_temp` (`LID`,`read_index`,`position`,`ID`) USING BTREE,
  ADD KEY `Dmso_contig_IDX` (`LID`,`position`) USING BTREE,
  ADD KEY `select_unmod_indx` (`LID`,`contig`,`read_index`,`position`);

--
-- AUTO_INCREMENT for dumped tables
--

--
-- AUTO_INCREMENT for table `Basecall`
--
ALTER TABLE `Basecall`
  MODIFY `ID` int(10) UNSIGNED NOT NULL AUTO_INCREMENT;

--
-- AUTO_INCREMENT for table `Basecall_Peaks`
--
ALTER TABLE `Basecall_Peaks`
  MODIFY `ID` int(10) UNSIGNED NOT NULL AUTO_INCREMENT;

--
-- AUTO_INCREMENT for table `Centroids`
--
ALTER TABLE `Centroids`
  MODIFY `ID` int(11) NOT NULL AUTO_INCREMENT;

--
-- AUTO_INCREMENT for table `Centroid_BPP`
--
ALTER TABLE `Centroid_BPP`
  MODIFY `ID` int(11) NOT NULL AUTO_INCREMENT;

--
-- AUTO_INCREMENT for table `Centroid_Distance`
--
ALTER TABLE `Centroid_Distance`
  MODIFY `ID` int(11) NOT NULL AUTO_INCREMENT;

--
-- AUTO_INCREMENT for table `Centroid_Probabilities`
--
ALTER TABLE `Centroid_Probabilities`
  MODIFY `ID` int(11) NOT NULL AUTO_INCREMENT;

--
-- AUTO_INCREMENT for table `Centroid_Regions`
--
ALTER TABLE `Centroid_Regions`
  MODIFY `ID` int(11) NOT NULL AUTO_INCREMENT;

--
-- AUTO_INCREMENT for table `Centroid_Secondary`
--
ALTER TABLE `Centroid_Secondary`
  MODIFY `ID` int(11) NOT NULL AUTO_INCREMENT;

--
-- AUTO_INCREMENT for table `Centroid_Secondary_Interactions`
--
ALTER TABLE `Centroid_Secondary_Interactions`
  MODIFY `ID` int(11) NOT NULL AUTO_INCREMENT;

--
-- AUTO_INCREMENT for table `Clusters`
--
ALTER TABLE `Clusters`
  MODIFY `ID` int(11) NOT NULL AUTO_INCREMENT;

--
-- AUTO_INCREMENT for table `Gmm`
--
ALTER TABLE `Gmm`
  MODIFY `ID` int(10) UNSIGNED NOT NULL AUTO_INCREMENT;

--
-- AUTO_INCREMENT for table `Library`
--
ALTER TABLE `Library`
  MODIFY `ID` int(11) NOT NULL AUTO_INCREMENT;

--
-- AUTO_INCREMENT for table `Lof`
--
ALTER TABLE `Lof`
  MODIFY `ID` int(10) UNSIGNED NOT NULL AUTO_INCREMENT;

--
-- AUTO_INCREMENT for table `Modified`
--
ALTER TABLE `Modified`
  MODIFY `ID` int(11) NOT NULL AUTO_INCREMENT;

--
-- AUTO_INCREMENT for table `Peaks`
--
ALTER TABLE `Peaks`
  MODIFY `ID` int(10) UNSIGNED NOT NULL AUTO_INCREMENT;

--
-- AUTO_INCREMENT for table `RAW`
--
ALTER TABLE `RAW`
  MODIFY `ID` int(11) NOT NULL AUTO_INCREMENT;

--
-- AUTO_INCREMENT for table `RAW_MOD`
--
ALTER TABLE `RAW_MOD`
  MODIFY `ID` int(11) NOT NULL AUTO_INCREMENT;

--
-- AUTO_INCREMENT for table `Reactivity_full`
--
ALTER TABLE `Reactivity_full`
  MODIFY `ID` int(11) NOT NULL AUTO_INCREMENT;

--
-- AUTO_INCREMENT for table `Read_depth`
--
ALTER TABLE `Read_depth`
  MODIFY `ID` int(11) NOT NULL AUTO_INCREMENT;

--
-- AUTO_INCREMENT for table `Read_depth_full`
--
ALTER TABLE `Read_depth_full`
  MODIFY `ID` int(11) NOT NULL AUTO_INCREMENT;

--
-- AUTO_INCREMENT for table `Structure`
--
ALTER TABLE `Structure`
  MODIFY `ID` int(10) UNSIGNED NOT NULL AUTO_INCREMENT;

--
-- AUTO_INCREMENT for table `Structure_BPP`
--
ALTER TABLE `Structure_BPP`
  MODIFY `ID` int(11) NOT NULL AUTO_INCREMENT;

--
-- AUTO_INCREMENT for table `Structure_Control`
--
ALTER TABLE `Structure_Control`
  MODIFY `ID` int(10) UNSIGNED NOT NULL AUTO_INCREMENT;

--
-- AUTO_INCREMENT for table `Structure_Library`
--
ALTER TABLE `Structure_Library`
  MODIFY `ID` int(11) NOT NULL AUTO_INCREMENT;

--
-- AUTO_INCREMENT for table `Structure_Metric`
--
ALTER TABLE `Structure_Metric`
  MODIFY `ID` int(11) NOT NULL AUTO_INCREMENT;

--
-- AUTO_INCREMENT for table `Structure_Metric_Secondary`
--
ALTER TABLE `Structure_Metric_Secondary`
  MODIFY `ID` int(11) NOT NULL AUTO_INCREMENT;

--
-- AUTO_INCREMENT for table `Structure_Probabilities`
--
ALTER TABLE `Structure_Probabilities`
  MODIFY `ID` int(11) NOT NULL AUTO_INCREMENT;

--
-- AUTO_INCREMENT for table `Structure_Secondary`
--
ALTER TABLE `Structure_Secondary`
  MODIFY `ID` int(11) NOT NULL AUTO_INCREMENT;

--
-- AUTO_INCREMENT for table `Structure_Secondary_Interactions`
--
ALTER TABLE `Structure_Secondary_Interactions`
  MODIFY `ID` int(11) NOT NULL AUTO_INCREMENT;

--
-- AUTO_INCREMENT for table `Unmodified`
--
ALTER TABLE `Unmodified`
  MODIFY `ID` int(11) NOT NULL AUTO_INCREMENT;
COMMIT;

/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
