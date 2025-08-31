<!-- NavigationTree.vue -->
<template>
  <ul class="list-unstyled mb-0">
    <li v-for="(chapter, index) in navTree" :key="index" class="mb-3">
      <!-- 一级：加粗，大一号字体 -->
      <div class="fw-bold fs-5 text-dark mb-2">{{ chapter.title }}</div>

      <ul class="list-unstyled mb-0">
        <li v-for="(section, sIndex) in chapter.children" :key="sIndex" class="mb-1">
          <!-- 二级：正常粗细，黑色，正常大小，无缩进 -->
          <div class="text-dark">{{ section.title }}</div>

          <!-- 三级：缩进 + 灰色 + 正常字体大小 -->
          <ul class="list-unstyled ms-3 mb-0 text-secondary">
            <li v-for="(sub, ssIndex) in section.children" :key="ssIndex">
              {{ sub.title }}
            </li>
          </ul>
        </li>
      </ul>
    </li>
  </ul>
</template>

<script>
export default {
  name: 'NavigationTree',
  data() {
    return {
      navTree: [
        {
          title: "引言",
          children: [
            {
              title: "背景介绍",
              children: [
                { title: "历史发展" },
                { title: "这是一个比较长的三级标题" },
                { title: "研究现状" },
                { title: "主要挑战" }
              ]
            },
            {
              title: "目标与结构",
              children: [
                { title: "本书目标" },
                { title: "章节安排" }
              ]
            }
          ]
        },
        {
          title: "核心概念",
          children: [
            {
              title: "基础理论",
              children: [
                { title: "定义与术语" },
                { title: "基本假设" }
              ]
            },
            {
              title: "关键技术",
              children: [
                { title: "算法原理" },
                { title: "实现方式" }
              ]
            },
            {
              title: "关键技术",
              children: [
                { title: "算法原理" },
                { title: "实现方式" }
              ]
            },
            {
              title: "关键技术",
              children: [
                { title: "算法原理" },
                { title: "实现方式" }
              ]
            },
            {
              title: "关键技术",
              children: [
                { title: "算法原理" },
                { title: "实现方式" }
              ]
            },
          ]
        }
      ]
    };
  }
};
</script>

<style scoped>
li {
  line-height: 1.5;
}
</style>
