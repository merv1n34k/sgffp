import { defineConfig } from 'vitepress'

export default defineConfig({
  title: 'sgffp',
  description: 'Read, write, and manipulate SnapGene files in Python',
  base: '/sgffp/',

  head: [
    ['link', { rel: 'icon', type: 'image/svg+xml', href: '/sgffp/logo.svg' }],
  ],

  themeConfig: {
    nav: [
      { text: 'Guide', link: '/guide/getting-started' },
      { text: 'Reference', link: '/reference/sgff-object' },
      { text: 'CLI', link: '/cli' },
      { text: 'Format Spec', link: '/format-spec' },
      {
        text: 'v0.15.0',
        items: [
          { text: 'Changelog', link: 'https://github.com/merv1n34k/sgffp/releases' },
          { text: 'PyPI', link: 'https://pypi.org/project/sgffp/' },
        ],
      },
    ],

    sidebar: {
      '/guide/': [
        {
          text: 'Guide',
          items: [
            { text: 'Getting Started', link: '/guide/getting-started' },
            { text: 'Core Concepts', link: '/guide/core-concepts' },
            { text: 'Reading & Writing', link: '/guide/reading-writing' },
            { text: 'History & Operations', link: '/guide/history-operations' },
          ],
        },
      ],
      '/reference/': [
        {
          text: 'API Reference',
          items: [
            { text: 'SgffObject & I/O', link: '/reference/sgff-object' },
            { text: 'Sequence & Features', link: '/reference/sequence-features' },
            { text: 'History', link: '/reference/history' },
            { text: 'Other Models', link: '/reference/other-models' },
            { text: 'Extras', link: '/reference/extras' },
            { text: 'Internals', link: '/reference/internals' },
          ],
        },
      ],
    },

    socialLinks: [
      { icon: 'github', link: 'https://github.com/merv1n34k/sgffp' },
    ],

    search: {
      provider: 'local',
    },

    editLink: {
      pattern: 'https://github.com/merv1n34k/sgffp/edit/master/docs/:path',
    },

    footer: {
      message: 'Released under the MIT License.',
    },
  },
})
